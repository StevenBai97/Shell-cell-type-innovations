# Run Cell Ranger count pipeline
cellranger count --id T --transcriptome ref/Cgi --fastqs ../scRNA/Cgigas/data/Trochophore --sample T --create-bam true --localcores 10
