# Identify reciprocal best hits (RBH) between human (Homo sapiens) and annelid worm (Pristina leidyi) proteins
mmseqs easy-rbh Homo_sapiens.GRCh38.pep.longest.fa ../Pleidyi/cellranger/ref/Ple.pep.fasta Hsa_Ple.rbh.tsv . --threads 20

# Extract worm orthologs of G2/M phase genes based on RBH results
grep -wf g2m_genes.id Hsa_Ple.rbh.tsv | cut -f 2 > g2m_genes.Ple

# Extract worm orthologs of S phase genes based on RBH results
grep -wf s_genes.id Hsa_Ple.rbh.tsv | cut -f 2 > s_genes.Ple
