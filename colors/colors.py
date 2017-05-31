colors = {}

for line in open('/home/unix/csmillie/box/tableau_colors.txt'):
    line = line.rstrip().split('\t')
    color = line[0]
    palette = line[1]
    order = line[2].split(': ')[1]
    if palette not in colors:
        colors[palette] = []
    colors[palette].append([order, color])

for palette in colors:
    colors[palette] = [xi[1] for xi in sorted(colors[palette])]

