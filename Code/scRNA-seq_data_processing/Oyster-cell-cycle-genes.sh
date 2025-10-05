# Identify reciprocal best hits (RBH) between human (Homo sapiens) and oyster (Crassostrea gigas) proteins
mmseqs easy-rbh Homo_sapiens.GRCh38.pep.longest.fa ../Cgigas/cellranger/ref/Cgi.pep.fasta Hsa_Cgi.rbh.tsv . --threads 20

# Extract oyster orthologs of G2/M phase genes based on RBH results
grep -wf g2m_genes.id Hsa_Cgi.rbh.tsv | cut -f 2 > g2m_genes.Cgi

# Map oyster G2/M orthologs to gene IDs using GFF3 annotation
grep -f g2m_genes.Cgi Cgi.gene.gff3 | sed 's/.*;gene=//g' | sed 's/;.*//g' | sort -u > g2m_genes.Cgi.geneid

# Extract oyster orthologs of S phase genes based on RBH results
grep -wf s_genes.id Hsa_Cgi.rbh.tsv | cut -f 2 > s_genes.Cgi

# Map oyster S phase orthologs to gene IDs using GFF3 annotation
grep -f s_genes.Cgi Cgi.gene.gff3 | sed 's/.*;gene=//g' | sed 's/;.*//g' | sort -u > s_genes.Cgi.geneid
