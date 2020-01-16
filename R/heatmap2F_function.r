#' heatmap2.F.R 
#'
#' This function runs differential gene analysis (linear modeling on normalized counts)
#' 
#' @keywords heatmap
#' @export
#' @import gplots
#' @import dynamicTreeCut
#' @import WGCNA
#' @examples
#' 
# AUTHOR: Freddie
# Note on naming:
# This function was originally referred to as 'heatmap.3',
# as an addon to 'heatmap.2'. However, a different 'heatmap.3'
# function has since been released, so we will now refer to this
# function as 'heatmap.F' for Freddy. ;)
#
# Instructions:
# The main purpose of this custom function is to make it easier 
# to change the color scale, to choose different distance
# metrics and clustering methods, which is awkward in existing
# heatmap functions. It also simplifies the process of choosing
# whether to cluster over rows, columns, both or neither.
#
# Arguments:
# data: The first argument is a numeric matrix containg the data to be
# plotted. The user should know beforehand whether the matrix
# contains both negative and positive values (e.g. logFC values), 
# or only positive (expression levels).
# 
# colors: This is a vector containing the colors to be plotted.
# The user can include any number of colors, but the middle one
# (e.g.  3 of 5, or 5 of 9) will correspond to values of 0 and the
# first colors to negative values.
# The color scale will then be adapted to the values in the data.
# If the data goes from e.g. 100 to -1, 100 will be represented by
# the strongest positive color, but -1 will look very close to 0.
# Defaults to 'darkred', 'orange', 'white', 'dodgerblue', 'darkblue'.
# This should show up well to the red-green color blind. 
#
# distmethod: This argument determines how distance is calculated.
# The user can pick any method compatible with heatmap.2, or one 
# of the following options:
# 'pearson' uses 1-pearson correlation as distance.
# 'spearman' uses 1-spearman correlation as distance.
# 'direct' uses the data matrix as a distance matrix directly.
# 'euclidean' uses eucludean distance.
# If anyone has ideas for further methods, they can easily be added.
# Defaults to 'pearson'.
#
# clustermethod: This argument determines how distances between
# genes are organized by hierarchical clustering. It can take any
# argument compatible with hclust. I think that heatmap.2 is only
# compatible with the default hclust method.
# Defaults to 'ward'
#
# clusterdim: Specifies whether to cluster data along rows, columns
# neither or both. Defaults to 'both'
#
# Version 2 has implemented rowside colors (cutoffmethod, cutoff).
# Version 3 will normalize the color scale to value distribution in data.

# library(gplots)
# library(dynamicTreeCut)
# library(WGCNA)

heatmap.F.4 = function(dataM, 
                       colors=c('darkblue', 'mediumblue', 'dodgerblue', 'white', 'orange', 'red', 'darkred'),
                       colorsaturation=0.25,
                       distmethod='pearson',
                       clustermethod='ward.D2',
                       clusterdim='row',
                       cutoffmethod='depth',
                       cutoff=3,
					   labRow=FALSE,
					   margins=c(20,2),
					   symbreaks=T,
					   scale=c("none"),
					   keytitle="log2 Fold Change",
                       main=NULL)
{
  if(length(colnames(dataM))==0){colnames(dataM) = as.character(1:ncol(dataM))}
  
  distmethod = match.arg(arg=distmethod, choices=c('bicor', 'pearson', 'spearman', 'direct', 'euclidean'))
  clusterdim = match.arg(arg=clusterdim, choices=c('both', 'row', 'column', 'none'))
  cutoffmethod = match.arg(arg=cutoffmethod, choices=c('depth', 'height', 'number'))
  
  if(clusterdim=='row'){Rowv=T; Colv=F; dendrogram='row'}
  if(clusterdim=='both'){Rowv=T; Colv=T; dendrogram='both'}
  if(clusterdim=='column'){Rolv=F; Colv=T; dendrogram='column'}
  if(clusterdim=='none'){Rowv=F; Colv=F; dendrogram='none'}
  
  # Color scale:	
  neg.col = rev(colors[1:ceiling(length(colors)/2)])
  pos.col = colors[ceiling(length(colors)/2):length(colors)]
  
  bins = diff(round(quantile(abs(dataM), seq(from=0, to=1, length.out=length(pos.col))^colorsaturation)/max(abs(range(dataM))), digits=3)*1000)
  a1 = list()
  a2 = list()
  for(i in 1:length(bins)){
    a1[[i]] = t(colorRamp(neg.col[c(i, i+1)])(seq(from=0, to=1, length.out=bins[i])))
    a2[[i]] = t(colorRamp(pos.col[c(i, i+1)])(seq(from=0, to=1, length.out=bins[i])))
  }
  a1 = matrix(unlist(a1), ncol=3, byrow=T); a2 = matrix(unlist(a2), ncol=3, byrow=T)
  
  b1 = rgb(red=a1[,1], green=a1[,2], blue=a1[,3], max=255)
  b2 = rgb(red=a2[,1], green=a2[,2], blue=a2[,3], max=255)
  
  #c = abs(range(dataM))
  #if(c[1]>c[2]){
  #	b2 = b2[1:round(length(b2)*c[2]/c[1])]
  #} else {
  #	b1 = b1[1:round(length(b1)*c[1]/c[2])]
  #}
  color.vector = c(rev(b1), b2)
  
  # Distance metric
  if(distmethod=='bicor'){fd = function(x){return(as.dist(1-bicor(t(x), use = 'pairwise.complete.obs')))}}
  if(distmethod=='pearson'){fd = function(x){return(as.dist(1-cor(t(x), method='pearson')))}}
  if(distmethod=='spearman'){fd = function(x){return(as.dist(1-cor(t(x), method='spearman')))}}
  if(distmethod=='direct'){fd = function(x){return(as.dist(x))}}
  if(distmethod=='euclidean'){fd = function(x){return(dist(x))}}
  
  # Clustering method
  fh = function(x){return(stats::hclust(x,method=clustermethod))}
  
  # Rowside colors
  
  if(cutoffmethod=='depth'){fc = function(M){return(cutreeHybrid(dendro=fh(fd(M)), distM=as.matrix(fd(M)), deepSplit=cutoff, verbose=0)$labels)}}
  if(cutoffmethod=='height'){fc = function(M){return(cutree(fh(fd(M)), h=cutoff))}}
  if(cutoffmethod=='number'){fc = function(M){return(cutree(fh(fd(M)), k=cutoff))}}
  
  if(dendrogram%in%c('none','column')){
    rowcol=rep('grey70', nrow(dataM))
  }else if(dendrogram%in%c('row','both')){
    rowcol=c('blue','red','orange','skyblue','yellow','black', 'darkblue','darkred','pink','purple','green','darkgreen','magenta','gray','brown','cyan','lightgreen','lightyellow','royalblue','white',"turquoise","hotpink")[suppressWarnings(fc(dataM))+1]
  }
  
  hm = heatmap.2(dataM,
                 col=color.vector,
                 hclustfun=fh,
                 distfun=fd,
                 key.title=NA,
                 keysize=1.5,
                 key.par=list(mar=c(10,1,1,1)),
                 key.xlab=keytitle,
                 trace='none',
                 density.info="none",
                 margins=margins,
                 labRow=labRow,
                 Rowv=Rowv,
                 Colv=Colv,
                 cexCol=1.5,
                 cexRow=0.5,
			         	 scale=scale,
                 dendrogram=dendrogram,
                 RowSideColors=rowcol,
                 symbreaks=symbreaks,
				 main=main)
  
  names(rowcol) = rownames(dataM)
  return(list(modules=rowcol[hm$rowInd], clustermatrix=t(hm$carpet)))
}
