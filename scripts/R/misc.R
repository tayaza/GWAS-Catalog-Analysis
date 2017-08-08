library(dendextend)
library(cluster)
library(factoextra)
library(NbClust)
library(pvclust)

#Measure distance of matrix
dist_eu <- dist(trait.matrix, method = "euclidean") 
fviz_dist_eu <- get_dist(trait.matrix, method = "euclidean", stand = FALSE) #snpdist_eu, snpdist_eu_st
fviz_dist(fviz_dist_eu, gradient = list(low = "red", mid = "white", high = "blue"))

#Assess clustering tendency

clustend <- get_clust_tendency(scale(trait.matrix), 100) #snp_clustend
clustend$hopkins_stat #[1] genes = 0.06561962; snps = 0.06702028; ctrl = 0.07315522 highly clusterable; well below threshold 0.5
clustend$plot #same as fviz_dist with stand=T

#Compute optimal number of clusters
#library(NbClust)
gap_sta <- clusGap(trait.matrix, FUN = hcut, K.max = 100, B = 100)
fviz_gap_stat(gap_sta)

#Compute hierarchical clustering
hc <- hclust(dist_eu, method = "complete") 
hc_w <- hclust(dist_eu, method = "ward.D2") #snp_hc_w, ctrl_hc_w
grp <- cutree(hc_w, k=22) #snp_grp, ctrl_grp
table(grp)

#Print .txt files containing traits in each cluster
for (i in 1:22){
  sink(paste("cluster", as.character(i), ".txt", sp=""), append=F, split=F)
  print(rownames(trait.matrix)[grp==i])
  sink()
}

#Plot grouped tree
png("grouped_tree.png", width = 10000, height = 10000, units = "px", pointsize = 22)
plot(hc_w, cex = 0.5)
rect.hclust(hc_w, k = 22, 
            border = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999"))
dev.off()


#Plot bootstrapped tree
set.seed(123) #Ensure reproducibility
res.pv <- pvclust(trait.matrix, method.hclust = "ward.D2", method.dist="euclidean", nboot = 1000) #snp_res.pv
png("tree_name.png", width = 10000, height = 10000, units = "px", pointsize = 22)
plot(res.pv1k, hang = -1, cex = 0.5)
pvrect(res.pv1k)
dev.off()



#corrplot(as.matrix(dist_eu), iscorr = FALSE, method = "color", order = "hclust", type = "upper")
#plot(hc_w)

#dend <- as.dendrogram(hc_w)
#dend <- rotate(dend)
#dend <- color_branches(dend)
#plot(dend)