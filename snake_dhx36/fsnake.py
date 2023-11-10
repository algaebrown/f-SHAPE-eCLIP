#configfile: "config.yaml"
# snakemake -j 30 -s fsnake.py --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn=1 -q home-yeo -e {params.error_out_file} -o /dev/null"
# snakemake -j 16 -s fsnake.py --cluster "qsub -l walltime={params.run_time}:00:00 -l nodes=1:ppn={params.cores} -q home-yeo" --directory=/home/hsher/rg4_seq/fshape_eclip_pipe/fshape_snake_pipe/

import pandas as pd
import os
import sys
import glob
from snake_config import CHROM_SIZES, GENOME_FA, GENOME_DICT, PYTHON3_PATH, pairs


include: "snake_config.py"

manifest = pd.read_table(MANIFEST, index_col = False, sep = ',')
sample_labels = manifest.Sample.tolist()



rule all:
    input:
        expand("reactivity/{comparison}..success.txt", 
        comparison=list(pairs.keys()))+
        expand('cov/{sample_label}.reverse.cov', sample_label = sample_labels)+
        expand('CIMS/{sample_label}.minus.bw', sample_label = sample_labels)
    output:
        "snakeCLIP.txt"
    params:
        error_out_file = "error_files/all",
        run_time = "1",
        cores = "1",
        memory = "20",
        job_name = "all"
    shell:
        "echo $(date)  > {output};"
        "echo created by Evan Boyle and the Yeo lab >> {output}"



rule header_free_sam:
    input:
        bam=lambda wildcards: glob.glob(manifest.loc[manifest.Sample == wildcards.sample_label]["bam"].values[0]),
    output:
        sam="processed_bam/{sample_label}.sam",
    params:
        run_time=6,
        error_out_file = "error_files/extract_read2",
        cores = "1",
    shell:
        """
        module load samtools;
        samtools view {input.bam} > {output.sam};
        """
rule split_reads:
    input:
        sam="processed_bam/{sample_label}.sam",
        
    output:
        mdtag="processed_bam/{sample_label}.mdtag",
        
    params:
        run_time=6,
        error_out_file = "error_files/pileup",
        cores = "1",
    shell:
        """
        # make strand specific
        python /home/hsher/f-SHAPE-eCLIP/snake_dhx36/splitRead.py {input.sam} {output.mdtag}
        
  
        """

rule sequece_dict:
    input:
        fa=GENOME_FA
    output:
        GENOME_DICT
    params:
        run_time=6,
        error_out_file = "error_files/pileup",
        fa=GENOME_FA,
        cores = "1",
    shell:
        """
        module load picard;
        picard CreateSequenceDictionary R={input.fa} O={output}
        """

rule split_reads_gatk:
    input:
        bam=lambda wildcards: glob.glob(manifest.loc[manifest.Sample == wildcards.sample_label]["bam"].values[0]),
        fa_dict=GENOME_DICT
    output:
        bam='processed_bam/{sample_label}.splicesplit.bam'
    params:
        run_time=6,
        error_out_file = "error_files/pileup",
        fa=GENOME_FA,
        cores = "3",
    shell:
        """
        # make strand specific
        module load gatk;
        gatk SplitNCigarReads \
            -R {params.fa} \
            -I {input.bam} \
            -O {output.bam}
  
        """

rule strand_specific_bam:
    input:
        bam='processed_bam/{sample_label}.splicesplit.bam'
    output:
        forward_bam="processed_bam/{sample_label}.splicesplit.forward.bam",
        reverse_bam="processed_bam/{sample_label}.splicesplit.reverse.bam"
    params:
        run_time=6,
        error_out_file = "error_files/pileup",
        cores = "2",
    shell:
        """
        # make strand specific
        module load samtools
        
        samtools view -h -F 0x10 {input.bam} | samtools view -Sb - > {output.forward_bam}
        samtools view -h -f 0x10 {input.bam}| samtools view -Sb - > {output.reverse_bam}
        """

rule sort_bam:
    input:
        forward_bam="processed_bam/{sample_label}.splicesplit.forward.bam",
        reverse_bam="processed_bam/{sample_label}.splicesplit.reverse.bam"
    output:
        forward_bam="processed_bam/{sample_label}.splicesplit.forward.sorted.bam",
        reverse_bam="processed_bam/{sample_label}.splicesplit.reverse.sorted.bam"
    params:
        run_time=6,
        error_out_file = "error_files/pileup",
        cores = "1",
    shell:
        """
        module load samtools
        samtools sort -o {output.forward_bam} {input.forward_bam}
        samtools sort -o {output.reverse_bam} {input.reverse_bam}
        """




rule pileup_mutation:
    input:
        forward_bam="processed_bam/{sample_label}.splicesplit.forward.sorted.bam",
        reverse_bam="processed_bam/{sample_label}.splicesplit.reverse.sorted.bam"
    output:
        forward_pileup='processed_bam/{sample_label}.forward.pileup',
        reverse_pileup='processed_bam/{sample_label}.reverse.pileup',

    params:
        run_time=16,
        error_out_file = "error_files/pileup",
        genomefa=GENOME_FA,
        cores = "6",
    shell:
        """
            module load samtools;
            samtools mpileup -s  {input.forward_bam} -f {params.genomefa} > {output.forward_pileup}
            samtools mpileup -s  {input.reverse_bam} -f {params.genomefa} > {output.reverse_pileup}
        """
    

rule CIMS_bw:
    input:
        plus_pileup='processed_bam/{sample_label}.forward.pileup',
        minus_pileup='processed_bam/{sample_label}.reverse.pileup',
    output:
        plus="CIMS/{sample_label}.plus.bw",
        minus="CIMS/{sample_label}.minus.bw",
        
    params:
        run_time=12,
        chr_size=CHROM_SIZES,
        python=PYTHON3_PATH,
        error_out_file = "error_files/CIMS_bw",
        cores = "1",

    shell:
        """
        # coverage
        {params.python} /home/hsher/projects/make_tracks/pileupToMismatchBw.py {input.plus_pileup} {params.chr_size} {output.plus}
        {params.python} /home/hsher/projects/make_tracks/pileupToMismatchBw.py {input.minus_pileup} {params.chr_size} {output.minus}
        """
rule pileup_to_cov:
    input:
        plus_pileup='processed_bam/{sample_label}.forward.pileup',
        minus_pileup='processed_bam/{sample_label}.reverse.pileup'
    output:
        plus='cov/{sample_label}.forward.cov',
        minus='cov/{sample_label}.reverse.cov'
    params:
        run_time=4,
        error_out_file = "error_files/CIMS_bw",
        cores = "1",
    shell:
        """
        cut -d $'\t' -f 1,2,3 {input.plus_pileup} > {output.plus};
        cut -d $'\t' -f 1,2,3 {input.minus_pileup} > {output.minus}
        """
        

rule sort_mdtag:
    input:
        mdtag="processed_bam/{sample_label}.mdtag",
    output:
        mdtag="processed_bam/{sample_label}.so.mdtag",
        
    params:
        run_time=6,
        error_out_file = "error_files/bam_pileup",
        cores="1"
    shell:
        """
        sort -k 1,1 -k 3,3n {input.mdtag} > {output.mdtag}
        """
rule countMutation: # TODO replace with pileup
    input:
        mdtag="processed_bam/{sample_label}.so.mdtag",
    output:
        cov="{sample_label}/none.neg.cov", 
        mut="{sample_label}/none.neg.mut", # TODO: very hard to see if it is successfully done.
       
    params:
        outdir="{sample_label}/",
        run_time=6,
        chr_size=CHROM_SIZES,
        python=PYTHON3_PATH,
        error_out_file = "error_files/CIMS_bw",
        cores="3"

    shell:
        """
        # coverage
        python /home/hsher/f-SHAPE-eCLIP/snake_dhx36/outputMut.py {input.mdtag} {params.outdir}
        """

rule prepare_bed:
    input:
        bed=BED,
    output:
        BED.replace('.bed', '.named.bed')
    
    params:
        run_time=2,
        python=PYTHON3_PATH,
        error_out_file = "error_files/namebed_bw",
        cores="1"

    shell:
        """
        # coverage
        {params.python} /home/hsher/f-SHAPE-eCLIP/snake_dhx36/name_bed.py {input.bed} {output}
        """


rule bed_reactivity:
    input:
        bed=BED.replace('.bed', '.named.bed'),
        treated=lambda wildcard: expand("{untreated_sample_label}/none.neg.cov",
            untreated_sample_label =  pairs[wildcard.comparison][0]),
        untreated=lambda wildcard: expand("{untreated_sample_label}/none.neg.cov",
            untreated_sample_label =  pairs[wildcard.comparison][1])
     # TODO: very hard to see if it is successfully done.
    output:
        outdir="reactivity/{comparison}..success.txt"
    params:
        run_time=6,
        fa=GENOME_FA,
        outdir="reactivity/{comparison}",
        python_path=PYTHON3_PATH,
        treated_dir=lambda wildcard: ','.join(pairs[wildcard.comparison][0]),
        untreated_dir=lambda wildcard: ','.join(pairs[wildcard.comparison][1]),
        cores="3"

    shell:
        """
        {params.python_path} /home/hsher/f-SHAPE-eCLIP/bedReactivities.py  -i {input.bed} \
            -g {params.fa} \
            -a {params.treated_dir} \
            -b {params.untreated_dir} \
            -o {params.outdir} 

        """