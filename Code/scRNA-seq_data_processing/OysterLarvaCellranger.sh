# Run Cell Ranger count pipeline
cellranger count --id T --transcriptome ref/Cgi --fastqs ../scRNA/Cgigas/data/Trochophore --sample T --create-bam true --localcores 10 &
cellranger count --id G4 --transcriptome ref/Cgi --fastqs ../scRNA/Cgigas/data/Gastrula --sample G4 --create-bam true --localcores 10 &
cellranger count --id G3 --transcriptome ref/Cgi --fastqs ../scRNA/Cgigas/data/Gastrula --sample G3 --create-bam true --localcores 10 &
cellranger count --id G2 --transcriptome ref/Cgi --fastqs ../scRNA/Cgigas/data/Gastrula --sample G2 --create-bam true --localcores 10 &
cellranger count --id G1 --transcriptome ref/Cgi --fastqs ../scRNA/Cgigas/data/Gastrula --sample G1 --create-bam true --localcores 10 &
