PYTHON3_PATH='~/miniconda3/bin/python'
MANIFEST ='/home/hsher/f-SHAPE-eCLIP/snake_dhx36/input/shape_samples.csv'

# Resources
CHROM_SIZES = "/home/hsher/gencode_coords/hg38.chrom.sizes"
GENOME_FA="/projects/ps-yeolab4/genomes/GRCh38/chromosomes/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
GENOME_DICT="/projects/ps-yeolab4/genomes/GRCh38/chromosomes/GRCh38_no_alt_analysis_set_GCA_000001405.15.dict" # use picard to make

BED='/home/hsher/rg4_seq/fshape_eclip_pipe/peaks/DHX36_D1.repr.bed'

pairs = {
    'IP_comapre':[['DHX36_N1_IP','DHX36_N2_IP'],
                    ['DHX36_D1_IP','DHX36_D2_IP']
    ],
    'IN_comapre':[['DHX36_N1_IN','DHX36_N2_IN'],
                ['DHX36_D1_IN','DHX36_D2_IN']
    ]
    }