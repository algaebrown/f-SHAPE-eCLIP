#configfile: "config.yaml"
# snakemake -j 30 -s fsnake.py --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn=1 -q home-yeo -e {params.error_out_file} -o /dev/null"
# snakemake -j 5 -s fsnake.py --cluster "qsub -l walltime={params.run_time}:00:00 -l nodes=1:ppn=1 -q home-yeo"

import pandas as pd
import os
import sys
import glob
from snake_config import CHROM_SIZES, GENOME_FA, GENOME_DICT, PYTHON3_PATH, PAIRING


include: "snake_config.py"

if not os.path.exists(MANIFEST): make_meta(MANIFEST)
manifest = pd.read_table(MANIFEST, index_col = False, sep = ',')

if not os.path.exists(PAIRING): make_meta(PAIRING)
pairs = pd.read_table(PAIRING, index_col = False, sep = ',')

sample_labels = manifest.Sample.tolist()
treated=pairs.treat.tolist()
untreated=pairs.untreat.tolist()

rule all:
    input:
        expand("{treated_sample_label}.{untreated_sample_label}.react/HISTMINUS.rx", 
        treated_sample_label = treated, # TODO: will combine pairs that does not make sense
        untreated_sample_label = untreated)+
        expand('processed_bam/{sample_label}.forward.pileup', sample_label = sample_labels)
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
    shell:
        """
        # make strand specific
        python splitRead.py {input.sam} {output.mdtag}
        
  
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
        run_time=6,
        error_out_file = "error_files/pileup",
        genomefa=GENOME_FA,
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
        minus="CIMS/{sample_label}.minus.bw"
    params:
        run_time=6,
        chr_size=CHROM_SIZES,
        python=PYTHON3_PATH,
        error_out_file = "error_files/CIMS_bw",

    shell:
        """
        # coverage
        {params.python} pileupToMismatchBw.py {input.plus_pileup} {params.chr_size} {output.plus}
        {params.python} pileupToMismatchBw.py {input.minus_pileup} {params.chr_size} {output.minus}
        """
rule pileup_to_cov:
    input:
        plus_pileup='processed_bam/{sample_label}.forward.pileup',
        minus_pileup='processed_bam/{sample_label}.reverse.pileup'
    output:
        plus='cov/{sample_label}.forward.cov',
        minus='cov/{sample_label}.reverse.cov'
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

    shell:
        """
        # coverage
        python outputMut.py {input.mdtag} {params.outdir}
        """

rule bed_reactivity:
    input:
        bed=BED,
        treated="{treated_sample_label}/none.neg.cov",
        untreated="{untreated_sample_label}/none.neg.cov" # TODO: very hard to see if it is successfully done.
    output:
        outdir="{treated_sample_label}.{untreated_sample_label}.success.txt"
    params:
        run_time=6,
        fa=GENOME_FA,
        treated_dir="{treated_sample_label}",
        untreated_dir="{untreated_sample_label}",
        outdir="{treated_sample_label}.{untreated_sample_label}",
        python_path=PYTHON3_PATH

    shell:
        """
        {params.python_path} /home/hsher/f-SHAPE-eCLIP/bedReactivities.py  -i {input.bed} \
            -g {params.fa} \
            -a {params.treated_dir} \
            -b {params.untreated_dir} \
            -o {params.outdir} 

        """