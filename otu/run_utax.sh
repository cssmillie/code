INFILE=$1
OUTFILE=$2

echo "usearch8 -utax $INFILE -db ~/db/utax/rdp_16s.udb -taxconfs ~/db/utax/rdp_16s_fl.tc -tt ~/db/utax/rdp_16s.tt -utaxout $OUTFILE"

