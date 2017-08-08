library("corrplot")
trait.data <- read.delim("genes_04.txt", header = T)
xtraits <- sort(unique(trait.data$trait_x))
ytraits <- sort(unique(trait.data$trait_y))

trait.matrix <- matrix(, nrow = length(ytraits), ncol = length(ytraits), dimnames = list(ytraits, ytraits))

for (i in rownames(trait.matrix)){
  for (j in colnames(trait.matrix)){
    i.index <- which(trait.data$trait_x == i)
    for (x in i.index){
      if (trait.data[x,2] == j){
        trait.matrix[i,j] <- trait.data[x,5]
      }
        
    }
  }
}
trait.matrix <- replace(trait.matrix, is.na(trait.matrix), 0)


col5 <-  colorRampPalette(c("pink", "green", "#0d2956", "#4d76c6","yellow", "red","#660404"))

corrplot(trait.matrix, tl.cex = 0.3, col = col5(100),
         cl.lim = c(0,1), addgrid.col = NA, tl.col = "black", is.corr = FALSE, order="hclust")


