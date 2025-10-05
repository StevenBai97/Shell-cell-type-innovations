# Identify reciprocal best hits (RBH) between human (Homo sapiens) and chaetognath (Paraspadella gotoi) proteins
mmseqs easy-rbh Homo_sapiens.GRCh38.pep.longest.fa ../Pgotoi/Pgo.pep.fasta Hsa_Pgo.rbh.tsv . --threads 20

# Extract chaetognath orthologs of G2/M phase genes based on RBH results
grep -wf g2m_genes.id Hsa_Pgo.rbh.tsv | cut -f 2 > g2m_genes.Pgo

# Map chaetognath G2/M orthologs to gene IDs using GFF3 annotation
grep -wf g2m_genes.Pgo ../Pgotoi/Pgot_genInfo_rr.txt | cut -f 5 > g2m_genes.Pgo.geneid

# Extract chaetognath orthologs of S phase genes based on RBH results
grep -wf s_genes.id Hsa_Pgo.rbh.tsv | cut -f 2 > s_genes.Pgo

# Map chaetognath S phase orthologs to gene IDs using GFF3 annotation
grep -wf s_genes.Pgo ../Pgotoi/Pgot_genInfo_rr.txt | cut -f 5 > s_genes.Pgo.geneid
