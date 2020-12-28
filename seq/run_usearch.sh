# run usearch_global
# input arguments
# 1 = query (fasta)
# 2 = db (fasta)
# 3 = id (fraction)
# 4 = out (text)
# example:
# run_usearch.sh query.fasta db.fasta 0.97 out.blast6
usearch -usearch_global $1 -db $2 -id $3 -blast6out $4
