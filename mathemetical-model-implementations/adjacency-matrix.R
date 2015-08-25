require(gdata)
require(stringr)

#setting the memory limit to 1GB
memory.limit(size=102400)

#loading all data in 2D matrix
rows <- readLines("../facebook/facebook_combined.txt")
egonet_table <- cbind(user_id = trim(substr(rows, 0, gregexpr(pattern ='\\s', rows))),  neighbour = trim(substr(rows, gregexpr(pattern ='\\s', rows), nchar(rows))))

#removing variable
rm(rows)

per_adj_mat <- 500

num_adj_mat <- nrow(egonet_table)/per_adj_mat

#adjacency matrix
adj_egonet <- matrix(rep(0), nrow = per_adj_mat, ncol = nrow(egonet_table), byrow = TRUE )
colnames(adj_egonet) <- seq(1:nrow(egonet_table))

write.table(adj_egonet, file = "adj-matrix.csv", sep = "|" , row.names=TRUE, col.names=FALSE)

for (k in 1:num_adj_mat ){
  
  rownames(as.numeric(adj_egonet)) <- seq(k,k+500)
  
  #iterating through all the rows n the table
  for (i in k:nrow(egonet_table)){
    adj_egonet[as.numeric(egonet_table[i,"user_id"])+1,as.numeric(egonet_table[i,"neighbour"])] <- 1
    if(i==k) break
  }
  
  #writing the adjacency matrix in a file
  write.table(adj_egonet, file = "adj-matrix.csv", sep = "|" , row.names=TRUE, col.names=FALSE , append=TRUE)
  
  #removing temp adjacency matrix
  rm(adj_egonet)
}
