for x in `ls *contam.fastq`; do
    y=${x%.*.*}
    for host in archaea bacteria fungi protozoa viral; do
        echo "python ~/box/fastq_to_fasta.py $x | blastall -p blastn -d ~/aviv/db/refseq/${host}/${host}.fna -m 8 -e 10 -o ./${host}/$y.$host.blastn"
    done
done
