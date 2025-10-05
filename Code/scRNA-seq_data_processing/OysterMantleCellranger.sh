# Convert raw FASTQ files into 10x-compatible format
dealFQ -r1 Unknown_BV836-001X0002_good_1.fq.gz -r2 Unknown_BV836-001X0002_good_2.fq.gz -o analysis -k bmk -t 30

# Build Cell Ranger reference
cellranger mkref --genome=Cgi --fasta=Cgi.genome.fasta --genes=Cgi.cellranger.gtf --nthreads 30

# Run Cell Ranger count pipeline
cellranger count --id M --transcriptome ref/Cgi --fastqs ../scRNA/Cgigas/data/Mantle --sample M --create-bam true --localcores 10
