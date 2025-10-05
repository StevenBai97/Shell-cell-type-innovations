# Identify reciprocal best hits (RBH) between human (Homo sapiens) and annelid worm (Pristina leidyi) proteins
mmseqs easy-rbh Homo_sapiens.GRCh38.pep.longest.fa ../Pleidyi/cellranger/ref/Ple.pep.fasta Hsa_Ple.rbh.tsv . --threads 20

# Extract worm orthologs of G2/M phase genes based on RBH results
grep -wf g2m_genes.id Hsa_Ple.rbh.tsv | cut -f 2 > g2m_genes.Ple

# Map oyster G2/M orthologs to gene IDs using GFF3 annotation
grep -f g2m_genes.Ple Ple.gene.gff3 | sed 's/.*;gene=//g' | sed 's/;.*//g' | sort -u > g2m_genes.Ple.geneid

# Extract oyster orthologs of S phase genes based on RBH results
grep -wf s_genes.id Hsa_Ple.rbh.tsv | cut -f 2 > s_genes.Ple

# Map oyster S phase orthologs to gene IDs using GFF3 annotation
grep -f s_genes.Ple Ple.gene.gff3 | sed 's/.*;gene=//g' | sed 's/;.*//g' | sort -u > s_genes.Ple.geneid
