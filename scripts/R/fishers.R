library("corrplot")
trait.data <- read.delim("gene_matrix.txt", header = T)
t_ordered <- read.delim("traits_new.txt", header = F)
trait_ordered <- t_ordered$V1

#Reverse order of matrix rownames (only if needed)
trait_reversed <- trait_ordered
for (i in seq(1,length(trait_ordered))){
  trait_reversed[i] <- trait_ordered[(length(trait_ordered) + 1 - i)]
}

#Instatiate empty matrix
trait_fischer.matrix <- matrix(, nrow = length(trait_ordered), ncol = length(trait_ordered), 
                       dimnames = list(trait_ordered, trait_ordered))

#Polpulate matrix
for (i in rownames(trait_fischer.matrix)){
  for (j in colnames(trait_fischer.matrix)){
    i.index <- which(trait.data$trait_x == i)
    for (x in i.index){
      if (trait.data[x,2] == j){
        trait_fischer.matrix[i,j] <- trait.data[x,6] #To Reverse p-values; 1-trait.data[x,6]
      }
      
    }
  }
}


#trait_fischer.matrix <- replace(trait_fischer.matrix, is.na(trait_fischer.matrix), 1)

for (i in rownames(trait_fischer.matrix)){
  for (j in colnames(trait_fischer.matrix)){
    if (i==j){
      trait_fischer.matrix[i,j] <- 0
    }
    
  }
}


col10 <-  colorRampPalette(c("green","#006d2c",  "#7f0000", "yellow", "#2166ac")) 
col11 <-  colorRampPalette(c("#fddbc7", "#ef8a62", "#ef8a62", "#b2182b", "#67a9cf", "#2166ac"))
cola <- colorRampPalette(c( "#023858","#045a8d","#0570b0","#3690c0","#74a9cf","#a6bddb","#d0d1e6","#ece7f2","#fff7fb"))
col12 <-  colorRampPalette(c("#e25904","#d81515", "#f7f7f7","#ffffff","#525252", "#252525","#000000"))
corrplot(trait_fischer.matrix, tl.cex = 0.3, col = col10(40), method = "color",
         cl.lim = c(0,1), addgrid.col = NA, tl.col = "black")


corrplot(trait_fischer.matrix, tl.cex = 0.1, col = col12(400), method = "color",
         cl.lim = c(0,1),  order="hclust", addgrid.col = NA, tl.col = "black")

#corrRect.hclust(trait.matrix, k = 300, col = "red", lwd = 1, method = "ward.D2")

