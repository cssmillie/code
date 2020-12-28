for x in `ls *contam.fastq`; do
    y=${x%.*.*}
    for host in archaea bacteria fungi protozoa viral; do
        echo "bowtie2 -x ~/aviv/db/bowtie2/${host}/${host} -U $x -S ./${host}/$y.$host.sam"
    done
done
