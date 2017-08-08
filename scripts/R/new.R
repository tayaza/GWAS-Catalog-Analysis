library("corrplot")
trait.data <- read.delim("genes/corrplot/gene_matrix.txt", header = T)
t_names <- read.delim("genes/corrplot/traits.txt", header = F)
trait_names <- t_names$V1
xtraits <- sort(unique(trait.data$trait_x))
ytraits <- sort(unique(trait.data$trait_y))

trait.matrix <- matrix(, nrow = length(trait_names), ncol = length(trait_names), 
                      dimnames = list(trait_names, trait_names))

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