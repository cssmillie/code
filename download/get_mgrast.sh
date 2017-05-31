PROJECT=$1
wget -m -nd -np -A *.fna.gz -X */*/*/processed ftp://ftp.metagenomics.anl.gov/projects/$PROJECT/
