# $1 = alignment (*.mafft)
# $2 = output name
raxmlHPC -m GTRGAMMA -s $1 -n $2 -p $RANDOM
