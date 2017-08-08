library("corrplot")
trait.data <- read.delim("traits04.txt", header = T)
xtraits <- sort(unique(trait.data$trait_x))
ytraits <- sort(unique(trait.data$trait_y))

trait.matrix <- matrix(, nrow = 301, ncol = 301, dimnames = list(xtraits, ytraits))

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

# for (i in rownames(mymatrix)){
#   for (j in colnames(mymatrix)){
#     if (i == j){
#       mymatrix[i,j] <- 1
#     } else
#         mymatrix[i,j] <- runif(1, 0, 1)
#       }
#       
# }
col1 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white",
                           "cyan", "#007FFF", "blue","#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))
col3 <- colorRampPalette(c("blue", "pink", "red"))
col4 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","#7FFF7F",
                           "cyan", "#007FFF", "blue","#00007F"))
col5 <-  colorRampPalette(c("black", "blue", "green", "white","lightpink", "red","brown"))

corrplot(trait.matrix, tl.cex = 0.3, method = "square", 
         cl.lim = c(0,1), addgrid.col = NA, tl.col = "gray", is.corr = FALSE, order="hclust")




#corrplot(trait.matrix, tl.cex = 0.3, na.label = "square", na.label.col = "white")
#trait.matrix <- replace(trait.matrix, is.na(trait.matrix), 0)