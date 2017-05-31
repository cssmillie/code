INFILE=$1
N_CPUS=$2
echo "ssu-prep -y $INFILE ssu $N_CPUS";
for((i=1; i<$N_CPUS; i++)); do echo "ssu-align ssu/$INFILE.$i ssu/ssu.$i"; done > ssu.1.sh;
echo "ssu-align --merge $N_CPUS ssu/$INFILE.$N_CPUS ssu/ssu.$N_CPUS" > ssu.2.sh;
echo "ssu-mask --stk2afa ssu"