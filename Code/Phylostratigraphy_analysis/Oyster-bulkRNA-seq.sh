# Quailty control of illumina reads with fastp (v 0.23.4)
## (Use SRR334223 as an example)
fastp -c -i SRR334223_1.fq.gz -I SRR334223_2.fq.gz -o SRR334223_1.clean.fq.gz -O SRR334223_2.clean.fq.gz

# Salmon (v.1.10.1) mapping and quantification
salmon index -t Cgi.transcript.fasta -i transcripts --keepDuplicates
salmon quant -i transcripts -l A -p 5 -1 SRR334223_1.clean.fq.gz -2 SRR334223_2.clean.fq.gz -o TC

# TMM normalization
mkdir quant
ls */quant.sf > quant/quant_files.fofn
perl abundance_estimates_to_matrix.pl --est_method salmon --quant_files quant/quant_files.fofn --name_sample_by_basedir --gene_trans_map quant/gene_trans_map --out_prefix quant/salmon
