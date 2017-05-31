import argparse, tableau_colors

# read input arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', help='input table')
parser.add_argument('-H', help='header', action='store_true', default=False)
parser.add_argument('-s', help='field separator (default = tab)', default='\t')
args = parser.parse_args()

# get colors
t20 = tableau_colors.colors['Tableau 20']
t20 = t20[0::2] + t20[1::2]

# get itol information
h = []
n = 0
for line in open(args.i):
    line = line.rstrip().split(args.s)
    # get header
    if not h:
        if args.H:
            h = line[1:]
            continue
        else:
            h = map(str, range(1,len(line)))
    # test for negatives
    for li in line[1:]:
        if float(li) < 0:
            n = 1

# duplicate header
if n ==  1:
    n = 2*len(h)
    h = ['%s_pos' %(hi) for hi in h] + ['%s_neg' %(hi) for hi in h]

# write itol header
print '''DATASET_PIECHART
SEPARATOR SPACE
DATASET_LABEL label1
COLOR #000000'''
print 'FIELD_COLORS ' + ' '.join(t20[:n])
print 'FIELD_LABELS ' + ' '.join(h)
print 'LEGEND_TITLE Legend'
print 'LEGEND_SHAPES ' + ' '.join(['2']*len(h))
print 'LEGEND_COLORS ' + ' '.join(t20[:len(h)])
print 'LEGEND_LABELS ' + ' '.join(h)
print 'DATA'

# write itol data
f = 0
for line in open(args.i):
    line = line.rstrip().split(args.s)
    
    # skip header
    if f == 0:
        f += 1
        continue
    
    # fix negatives if necessary
    new = [0]*len(h)
    for i in range(1, len(line)):
        if float(line[i]) >= 0:
            new[i-1] = float(line[i])
        else:
            new[i + len(line) - 2] = abs(float(line[i]))
    
    # insert position and radius
    radius = sum(new)
    print ' '.join(map(str, [line[0], 0, radius] + new))
