require(RColorBrewer)

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

set.colors = c(brewer.pal(9, 'Set1'), brewer.pal(7, 'Set2'), brewer.pal(12, 'Set3')[c(3,9,8)], 'violetred4')
set.colors[6] = 'khaki2'
set.colors[8] = 'lightskyblue2'
set.colors = rep(set.colors, 10)

tsne.colors = c(brewer.pal(7, 'Set2'), brewer.pal(9, 'Set1'), brewer.pal(12, 'Set3')[c(3,9,8)], 'violetred4')
tsne.colors[13] = 'khaki2'
tsne.colors[15] = 'lightskyblue2'
tsne.colors = rep(tsne.colors, 10)


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
