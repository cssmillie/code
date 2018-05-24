require(RColorBrewer)

desat = function(cols, sat=0.5) {
    X = diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
    hsv(X[1,], X[2,], X[3,])
}

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

set.colors = c(brewer.pal(9, 'Set1'), brewer.pal(7, 'Set2'), brewer.pal(12, 'Set3')[c(3,9,8)], 'violetred4')
set.colors[6] = 'khaki2'
set.colors[8] = 'lightskyblue2'
set.colors = rep(set.colors, 10)

tsne.colors = c(brewer.pal(7, 'Set2'), brewer.pal(9, 'Set1'), brewer.pal(12, 'Set3')[c(3,9,8)], 'violetred4')
tsne.colors[13] = 'khaki2'
tsne.colors[15] = 'lightskyblue2'
tsne.colors = rep(tsne.colors, 10)
tsne.colors = desat(tsne.colors, sat=.75)

material.heat <- function(n)
{
    mh = c(
        #"#607D8B", #blue grey
        "#283593", #indigo 800
        "#3F51B5",  #indigo
        "#2196F3", # blue
        #"#03A9F4", # light blue
        "#00BCD4", #cyan
        #"#009688", # teal
        "#4CAF50", #green
        "#8BC34A", #light green
        "#CDDC39", # lime
        "#FFEB3B", #yellow
        "#FFC107", # amber
        "#FF9800", # organe
        "#FF5722", # deep orange)
        "#f44336")
    colorRampPalette(mh)(n)
}

nmf.colors = c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")

tol1qualitative=c("#4477AA")
tol2qualitative=c("#4477AA", "#CC6677")
tol3qualitative=c("#4477AA", "#DDCC77", "#CC6677")
tol4qualitative=c("#4477AA", "#117733", "#DDCC77", "#CC6677")
tol5qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677")
tol6qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499")
tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")
tol8qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677","#AA4499")
tol9qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")
tol10qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
tol11qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
tol12qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")


# material colors
# ---------------

material = read.table('~/code/single_cell/material_colors.txt', row.names=1, comment.char='', stringsAsFactors=F)
rownames(material) = tolower(rownames(material))

shades = function(color, n){
    colorRampPalette(material[color,])(n)
}
disc.colors = function(n, shade=6){
    disc.order = c('red', 'light_blue', 'amber', 'green', 'grey')
    if(n <= length(disc.order)){
        material[disc.order[1:n], shade]
    } else if(n <= nrow(material)){
        material[as.integer(seq(from=1, to=nrow(material), length.out=n)), shade]
    } else {
        colorRampPalette(material[,shade])(n)
    }
}
cont.colors = function(shade=6){
    colorRampPalette(material[,shade])(100)
}
name.colors = function(names, shade=6){
    material[names, shade]
}
