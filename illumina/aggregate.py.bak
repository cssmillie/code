# Aggregate output
print 'perl /seq/regev_genome_portal/SOFTWARE/KCO/RNASEQ_pipeline/aggregate_links_to_sample_outputs.pl ./'

# Run QC
print 'perl /seq/regev_genome_portal/SOFTWARE/KCO/RNASEQ_pipeline/util/summarize_rnaseqQC_results.pl samples.txt ./ RNASEQC_STAR > ./rnaseqQC.txt'

# Combine RSEM
print 'find ./RSEM_*/ -type f  | egrep "genes.results" > rsem.genes.list'
print 'perl /seq/regev_genome_portal/SOFTWARE/KCO/RNASEQ_pipeline/merge_RSEM_output_to_matrix.pl --rsem_files rsem.genes.list --mode counts > rsem.genes.counts.matrix'

