# Functional annotation using eggNOG-mapper
emapper.py -m diamond --tax_scope 33208 -i Cgi.pep.fasta --evalue 1e-5 --output Cgi.pep.emapper_out --data_dir db/eggnog-mapper-data/ --cpu 30
Rscript emapperx.R Cgi.pep.emapper_out.emapper.annotations Cgi.pep.fasta
