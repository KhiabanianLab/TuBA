BiclusteringFunction <- function(FileNameGenePairs,FileNameBinaryMatrix,OverlapCutOff)
{
  library(data.table)
  
  library(plyr)
  
  library(igraph)
  
  SignificanceLevel <- 10^(-OverlapCutOff)
  
  Gene.Pairs.df <- fread(FileNameGenePairs)
  
  Genes.Samples.Binary.df <- fread(FileNameBinaryMatrix)
  
  Non.Zero.Genes <- as.character(Genes.Samples.Binary.df$Gene.ID)
  
  Binary.Matrix.For.Genes.Outliers <- as.matrix(Genes.Samples.Binary.df[,-1])
  
  Sample.IDs <- colnames(Genes.Samples.Binary.df[,-1])
  
  rm(Genes.Samples.Binary.df)
  
  Gene.Pairs.Above.Significance.Level <- which(Gene.Pairs.df$p.val <= SignificanceLevel)
  
  Revised.Col1 <- Gene.Pairs.df$Gene.1[Gene.Pairs.Above.Significance.Level]
  
  Revised.Col2 <- Gene.Pairs.df$Gene.2[Gene.Pairs.Above.Significance.Level]
  
  All.Genes <- union(Gene.Pairs.df$Gene.1[Gene.Pairs.Above.Significance.Level],Gene.Pairs.df$Gene.2[Gene.Pairs.Above.Significance.Level])

  #Binary matrix that associates only those genes that are part of triangular cliques at the very least
  Pruning.Matrix <- matrix(0,nrow = length(Non.Zero.Genes),ncol = length(Non.Zero.Genes))
  for(i in 1:length(All.Genes))
  {
    TempI <- All.Genes[i]
    
    List.of.Associated.Genes1 <- Revised.Col2[which(Revised.Col1 == TempI)]
  
    List.of.Associated.Genes2 <- Revised.Col1[which(Revised.Col2 == TempI)]
  
    List.of.Associated.Genes <- unique(c(List.of.Associated.Genes1,List.of.Associated.Genes2))
  
    for (j in 1:length(List.of.Associated.Genes))
    {
    
      List.of.Associated.Genes1 <- Revised.Col2[which(Revised.Col1 == List.of.Associated.Genes[j])]
    
      List.of.Associated.Genes2 <- Revised.Col1[which(Revised.Col2 == List.of.Associated.Genes[j])]
    
      List.of.Associated.Genes.Temp <- unique(c(List.of.Associated.Genes1,List.of.Associated.Genes2))
    
      Shared.List.of.Associated.Genes <- intersect(List.of.Associated.Genes.Temp,List.of.Associated.Genes)
    
      Pruning.Matrix[Shared.List.of.Associated.Genes,TempI] <- 1
    }
  }

  Pruning.Matrix[lower.tri(Pruning.Matrix)] <- 0

  Pruned.Graph.Array <- which(Pruning.Matrix == 1,arr.ind = T)

  Pruned.Graph.Col1 <- Pruned.Graph.Array[,1]

  Pruned.Graph.Col2 <- Pruned.Graph.Array[,2]
  
  rm(Pruning.Matrix)

  Pruned.Graph.df <- data.frame(Pruned.Graph.Col1,Pruned.Graph.Col2)
  colnames(Pruned.Graph.df) <- c("Node1","Node2")

  Pruned.Graph.Unique.Genes <- unique(c(Pruned.Graph.df$Node1,Pruned.Graph.df$Node2))

  Degrees.Pruned.Graph.Unique.Genes <- vector(mode = "numeric",length = length(Pruned.Graph.Unique.Genes))
  for (i in 1:length(Pruned.Graph.Unique.Genes))
  {
    Temp.Length1 <- length(which(Pruned.Graph.df$Node1 == Pruned.Graph.Unique.Genes[i]))
    Temp.Length2 <- length(which(Pruned.Graph.df$Node2 == Pruned.Graph.Unique.Genes[i]))
    Degrees.Pruned.Graph.Unique.Genes[i] <- Temp.Length1 + Temp.Length2
  }

  Pruned.Graph.Unique.Genes <- Pruned.Graph.Unique.Genes[order(Degrees.Pruned.Graph.Unique.Genes,decreasing = T)]

  #Code to identify all distinct clusters in graph

  k=1
  List.of.Genes.In.Clusters <- vector(mode = "list")
  while (length(Pruned.Graph.Unique.Genes) >= 1)
  {
    GeneName <- Pruned.Graph.Unique.Genes[1]
  
    List.of.Associated.Genes1 <- Pruned.Graph.df$Node2[which(Pruned.Graph.df$Node1 == GeneName)]
  
    List.of.Associated.Genes2 <- Pruned.Graph.df$Node1[which(Pruned.Graph.df$Node2 == GeneName)]
  
    List.of.Associated.Genes <- union(List.of.Associated.Genes1,List.of.Associated.Genes2)
  
    Length.List.of.Associated.Genes <- length(List.of.Associated.Genes)
  
    Length.List.of.Associated.Genes1 <- length(List.of.Associated.Genes[-c(1,2)])
  
    while (Length.List.of.Associated.Genes/Length.List.of.Associated.Genes1 > 1)
    {
      Length.List.of.Associated.Genes1 <- length(List.of.Associated.Genes)
      for (i in 1:length(List.of.Associated.Genes))
      {
        List.of.Associated.Genes1 <- Pruned.Graph.df$Node2[which(Pruned.Graph.df$Node1 == List.of.Associated.Genes[i])]
      
        List.of.Associated.Genes2 <- Pruned.Graph.df$Node1[which(Pruned.Graph.df$Node2 == List.of.Associated.Genes[i])]
      
        List.of.Associated.Genes.Temp <- unique(c(List.of.Associated.Genes1,List.of.Associated.Genes2))
      
        List.of.Associated.Genes <- unique(c(List.of.Associated.Genes,List.of.Associated.Genes.Temp))
      
      }
      Length.List.of.Associated.Genes <- length(List.of.Associated.Genes)
    }
    List.of.Genes.In.Clusters[[k]] <- as.numeric(List.of.Associated.Genes)
    Pruned.Graph.Unique.Genes <- Pruned.Graph.Unique.Genes[!Pruned.Graph.Unique.Genes %in% List.of.Associated.Genes]
    k <- k+1
  }


  #Code to find seeds in the clusters

  List.of.Seeds <- vector(mode = "list")
  m <- 1
  t <- 1
  i <- 0
  Temp.Largest.Clique.Length <- 1
  x2 <- data.frame(Pruned.Graph.df$Node1[Pruned.Graph.df$Node1 %in% List.of.Genes.In.Clusters[[1]] & Pruned.Graph.df$Node2 %in% List.of.Genes.In.Clusters[[1]]],Pruned.Graph.df$Node2[Pruned.Graph.df$Node1 %in% List.of.Genes.In.Clusters[[1]] & Pruned.Graph.df$Node2 %in% List.of.Genes.In.Clusters[[1]]])
  colnames(x2) <- c("Node1","Node2")
  while (i < length(List.of.Genes.In.Clusters))
  {
    if (Temp.Largest.Clique.Length < 3 | length(x2$Node1) ==0)
    {
      i <- i+1
      x2 <- data.frame(Pruned.Graph.df$Node1[Pruned.Graph.df$Node1 %in% List.of.Genes.In.Clusters[[i]] & Pruned.Graph.df$Node2 %in% List.of.Genes.In.Clusters[[i]]],Pruned.Graph.df$Node2[Pruned.Graph.df$Node1 %in% List.of.Genes.In.Clusters[[i]] & Pruned.Graph.df$Node2 %in% List.of.Genes.In.Clusters[[i]]])
      colnames(x2) <- c("Node1","Node2")
    }
  
    g=graph.data.frame(x2,directed=FALSE)
  
    List.of.Largest.Cliques <- vector(mode = "list")
  
    List.of.Largest.Cliques <- largest_cliques(g)
  
    List.of.Seeds[[t]] <- as.numeric(V(g)$name[unlist(List.of.Largest.Cliques[[1]])])
  
    Nodes.To.Be.Removed <- List.of.Seeds[[t]]
  
    Temp.Sub.Cluster <- vector(mode = "numeric")
    List.Temp.Sub.Cluster <- vector(mode = "list")
    Temp.Intersection.Cliques <- vector(mode = "numeric")
  
    for (m in 1:length(List.of.Largest.Cliques))
    {
      Temp.Sub.Cluster <- as.numeric(V(g)$name[unlist(List.of.Largest.Cliques[[m]])])
      List.Temp.Sub.Cluster[[m]] <- Temp.Sub.Cluster
      Temp.Intersection.Cliques[m] <- length(intersect(Nodes.To.Be.Removed,Temp.Sub.Cluster))
    }
  
    List.of.Seeds[[t]] <- unique(c(List.of.Seeds[[t]],unlist(List.Temp.Sub.Cluster[which(Temp.Intersection.Cliques !=0)])))
  
    Nodes.To.Be.Removed <- unlist(List.of.Seeds)
  
    x2 <- data.frame(x2$Node1[!x2$Node1 %in% Nodes.To.Be.Removed & !x2$Node2 %in% Nodes.To.Be.Removed],x2$Node2[!x2$Node1 %in% Nodes.To.Be.Removed & !x2$Node2 %in% Nodes.To.Be.Removed])
    colnames(x2) <- c("Node1","Node2")
  
    t <- t+1
    if (length(x2$Node1) > 0)
    {
      g1=graph.data.frame(x2,directed=FALSE)
      Temp.Largest.Clique.Length <- length(as.numeric(V(g1)$name[unlist(largest_cliques(g1)[[1]])]))
    }
  }

  #Edgelist for seeds
  y1 <- c(0)
  y2 <- c(0)
  Lengths.Edgelist.y <- vector(mode = "numeric")
  for (i in 1:length(List.of.Seeds))
  {
    y1 <- c(y1,Pruned.Graph.df$Node1[Pruned.Graph.df$Node1 %in% List.of.Seeds[[i]] & Pruned.Graph.df$Node2 %in% List.of.Seeds[[i]]])
    y2 <- c(y2,Pruned.Graph.df$Node2[Pruned.Graph.df$Node1 %in% List.of.Seeds[[i]] & Pruned.Graph.df$Node2 %in% List.of.Seeds[[i]]])
    Lengths.Edgelist.y[i] <- length(Pruned.Graph.df$Node1[Pruned.Graph.df$Node1 %in% List.of.Seeds[[i]] & Pruned.Graph.df$Node2 %in% List.of.Seeds[[i]]])
  }

  y <- data.frame(y1[-1],y2[-1])
  colnames(y) <- c("Node1","Node2")

  y$Node1 <- Non.Zero.Genes[y$Node1]
  y$Node2 <- Non.Zero.Genes[y$Node2]

  Vec.Edgelist.Seed.Labels <- vector(mode = "numeric")
  for (i in 1:length(Lengths.Edgelist.y))
  {
    Vec.Edgelist.Seed.Labels <- c(Vec.Edgelist.Seed.Labels,rep(i,Lengths.Edgelist.y[i]))
  }

  Seed.Edgelist.df <- data.frame(y$Node1,y$Node2,Vec.Edgelist.Seed.Labels)
  colnames(Seed.Edgelist.df) <- c("Node1","Node2","Seed.No")
  
  File.Name <- as.character(paste(gsub("SignificantGenePairs.csv","",FileNameGenePairs),"OverlapCutOff",toupper(as.character(SignificanceLevel)),"_Edgelist_Seeds.csv",sep = ""))
  
  write.table(Seed.Edgelist.df,file = File.Name,row.names = F,col.names = T,sep = ",")

  #Code to generate Seed versus Samples Matrix
  List.of.Samples.Per.Seed <- vector(mode = "list",length = length(List.of.Seeds))
  for (i in 1:length(List.of.Seeds))
  {
    Temp.List.of.Genes.in.Seed <- List.of.Seeds[[i]]
    x3 <- data.frame(Pruned.Graph.df$Node1[Pruned.Graph.df$Node1 %in% List.of.Seeds[[i]] & Pruned.Graph.df$Node2 %in% List.of.Seeds[[i]]],Pruned.Graph.df$Node2[Pruned.Graph.df$Node1 %in% List.of.Seeds[[i]] & Pruned.Graph.df$Node2 %in% List.of.Seeds[[i]]])
    colnames(x3) <- c("Node1","Node2")
    Temp.Gene.Table.Col1 <- x3$Node1
    Temp.Gene.Table.Col2 <- x3$Node2
  
    Temp.List.of.Samples.Per.Seed <- vector(mode = "list")
    for (j in 1:length(Temp.Gene.Table.Col1))
    {
      TempJ.Col1 <- Temp.Gene.Table.Col1[j]
      TempJ.Col2 <- Temp.Gene.Table.Col2[j]
      Temp.List.of.Samples.Per.Seed[[j]] <- intersect(which(Binary.Matrix.For.Genes.Outliers[TempJ.Col1,] == 1),which(Binary.Matrix.For.Genes.Outliers[TempJ.Col2,] == 1))
    }
    List.of.Samples.Per.Seed[[i]] <- unique(unlist(Temp.List.of.Samples.Per.Seed))
  }
  
  Seed.Sample.Signature.Matrix <- matrix(0,nrow = length(List.of.Seeds),ncol = length(Sample.IDs))
  
  for (i in 1:length(List.of.Seeds))
  {
    Seed.Sample.Signature.Matrix[i,List.of.Samples.Per.Seed[[i]]] <- 1
  }
  
  Seed.Nos <- seq(1,nrow(Seed.Sample.Signature.Matrix),1)
  
  Seed.Sample.Signature.df <- data.frame(Seed.Nos,Seed.Sample.Signature.Matrix)
  
  colnames(Seed.Sample.Signature.df) <- c("Seed.No",Sample.IDs)
  
  File.Name <- as.character(paste(gsub("SignificantGenePairs.csv","",FileNameGenePairs),"OverlapCutOff",toupper(as.character(SignificanceLevel)),"_SeedsSamplesMatrix.csv",sep = ""))
  
  write.table(Seed.Sample.Signature.df, file = File.Name,row.names = F,col.names = T,sep = ",")
  
  #Code to obtain names of genes per seed
  
  Seed.No.Vec <- vector(mode = "numeric")
  Genes.List.Seeds <- vector(mode = "character")
  for (i in 1:length(List.of.Seeds))
  {
    Seed.No.Temp.Vec <- rep(i,length(List.of.Seeds[[i]]))
    Seed.No.Vec <- c(Seed.No.Vec,Seed.No.Temp.Vec)
    Genes.Temp.List.Seeds <- as.character(Non.Zero.Genes[List.of.Seeds[[i]]])
    Genes.List.Seeds <- c(Genes.List.Seeds,Genes.Temp.List.Seeds)
  }

  Genes.List.Seeds.df <- data.frame(Genes.List.Seeds,Seed.No.Vec)
  colnames(Genes.List.Seeds.df) <- c("Gene.ID","Seed.No")
  
  File.Name <- as.character(paste(gsub("SignificantGenePairs.csv","",FileNameGenePairs),"OverlapCutOff",toupper(as.character(SignificanceLevel)),"_GeneNames_Seeds.csv",sep = ""))
  
  write.table(Genes.List.Seeds.df,file = File.Name,row.names = F,col.names = T,sep = ",")
  
  #Code for sowing the "seeds" in pruned graph

  x4 <- data.frame(c(0),c(0))
  colnames(x4) <- c("Node1","Node2")
  Lengths.Edgelist.x3 <- vector(mode = "numeric")
  New.List.of.Seeds <- vector(mode = "list",length = length(List.of.Seeds))
  for (i in 1:length(List.of.Seeds))
  {
    x3 <- data.frame(Pruned.Graph.df$Node1[Pruned.Graph.df$Node1 %in% List.of.Seeds[[i]] | Pruned.Graph.df$Node2 %in% List.of.Seeds[[i]]],Pruned.Graph.df$Node2[Pruned.Graph.df$Node1 %in% List.of.Seeds[[i]] | Pruned.Graph.df$Node2 %in% List.of.Seeds[[i]]])
    colnames(x3) <- c("Node1","Node2")
  
    Common.List.Vec <- c(x3$Node1,x3$Node2)
    while (length(which(count(Common.List.Vec)$freq <2)) > 0)
    {
      Common.List.Vec <- c(x3$Node1,x3$Node2)
    
      Vec.of.Leaves <- count(Common.List.Vec)$x[which(count(Common.List.Vec)$freq < 2)]
    
      x3 <- data.frame(x3$Node1[!x3$Node1 %in% Vec.of.Leaves & !x3$Node2 %in% Vec.of.Leaves],x3$Node2[!x3$Node1 %in% Vec.of.Leaves & !x3$Node2 %in% Vec.of.Leaves])
      colnames(x3) <- c("Node1","Node2")
    
      Common.List.Vec <- c(x3$Node1,x3$Node2)
    
    }
    Lengths.Edgelist.x3[i] <- length(x3$Node1)
    New.List.of.Seeds[[i]] <- unique(c(x3$Node1,x3$Node2))
    x4 <- data.frame(c(x4$Node1,x3$Node1),c(x4$Node2,x3$Node2))
    colnames(x4) <- c("Node1","Node2")
  }

  x4 <- x4[-1,]

  #List of Samples per Bicluster

  List.of.Samples.Per.Bicluster <- vector(mode = "list",length = length(New.List.of.Seeds))
  List.of.Genes.Per.Bicluster <- vector(mode = "list",length = length(New.List.of.Seeds))
  Degrees.of.Genes.In.Biclusters<- vector(mode = "list",length = length(New.List.of.Seeds))
  Bicluster.Nos.List <- vector(mode = "list",length = length(New.List.of.Seeds))
  k <- 0
  TempI <- 0
  for (i in 1:length(New.List.of.Seeds))
  {
    k <- TempI+1
    TempI <- TempI + Lengths.Edgelist.x3[i]
    Temp.Gene.Table.Col1 <- x4$Node1[k:TempI]
    Temp.Gene.Table.Col2 <- x4$Node2[k:TempI]
    List.of.Genes.Per.Bicluster[[i]] <- count(c(Temp.Gene.Table.Col1,Temp.Gene.Table.Col2))$x
    Degrees.of.Genes.In.Biclusters[[i]] <- count(c(Temp.Gene.Table.Col1,Temp.Gene.Table.Col2))$freq
    Bicluster.Nos.List[[i]] <- rep(i,length(List.of.Genes.Per.Bicluster[[i]]))
    Temp.List.of.Samples.Per.Bicluster <- vector(mode = "list")
    for (j in 1:Lengths.Edgelist.x3[i])
    {
      TempJ.Col1 <- Temp.Gene.Table.Col1[j]
      TempJ.Col2 <- Temp.Gene.Table.Col2[j]
      Temp.List.of.Samples.Per.Bicluster[[j]] <- intersect(which(Binary.Matrix.For.Genes.Outliers[TempJ.Col1,] == 1),which(Binary.Matrix.For.Genes.Outliers[TempJ.Col2,] == 1))
    }
    List.of.Samples.Per.Bicluster[[i]] <- unique(unlist(Temp.List.of.Samples.Per.Bicluster))
  }

  #Matrix with Biclusters and Samples

  Bicluster.Sample.Signature.Matrix <- matrix(0,nrow = length(New.List.of.Seeds),ncol = length(Sample.IDs))

  for (i in 1:length(New.List.of.Seeds))
  {
    Bicluster.Sample.Signature.Matrix[i,List.of.Samples.Per.Bicluster[[i]]] <- 1
  }

  Bicluster.Sample.Signature.df <- data.frame(seq(1,length(New.List.of.Seeds),1),Bicluster.Sample.Signature.Matrix)
  colnames(Bicluster.Sample.Signature.df) <- c("Bicluster.No",Sample.IDs)
  
  File.Name <- as.character(paste(gsub("SignificantGenePairs.csv","",FileNameGenePairs),"OverlapCutOff",toupper(as.character(SignificanceLevel)),"_BiclusterSamplesMatrix.csv",sep = ""))
  
  write.table(Bicluster.Sample.Signature.df,file = File.Name,row.names = F,col.names = T,sep = ",")
  
  #Table with genes and their respective degrees in each bicluser
  
  Genes.Degrees.Bicluster.df <- data.frame(Non.Zero.Genes[unlist(List.of.Genes.Per.Bicluster)],unlist(Degrees.of.Genes.In.Biclusters),unlist(Bicluster.Nos.List))
  colnames(Genes.Degrees.Bicluster.df) <- c("Gene.ID","Degree.In.Bicluster","Bicluster.No")
  
  File.Name <- as.character(paste(gsub("SignificantGenePairs.csv","",FileNameGenePairs),"OverlapCutOff",toupper(as.character(SignificanceLevel)),"_GenesDegreesInBiclusters.csv",sep = ""))
  
  write.table(Genes.Degrees.Bicluster.df,file = File.Name,row.names = F,col.names = T,sep = ",")
  
  #Edgelist
  Vec.Edgelist.Bicluster.Labels <- vector(mode = "numeric")
  for (i in 1:length(Lengths.Edgelist.x3))
  {
    Vec.Edgelist.Bicluster.Labels <- c(Vec.Edgelist.Bicluster.Labels,rep(i,Lengths.Edgelist.x3[i]))
  }
  
  Edgelist.GeneNos.df <- data.frame(x4$Node1,x4$Node2,Vec.Edgelist.Bicluster.Labels)
  colnames(Edgelist.GeneNos.df) <- c("Node1","Node2","Bicluster.No")

  Temp.Vec1 <- Non.Zero.Genes[Edgelist.GeneNos.df$Node1]
  Temp.Vec2 <- Non.Zero.Genes[Edgelist.GeneNos.df$Node2]

  Edgelist.df <- data.frame(Temp.Vec1,Temp.Vec2,Edgelist.GeneNos.df$Bicluster.No)
  colnames(Edgelist.df) <- c("Node1","Node2","Bicluster.No")
  
  File.Name <- as.character(paste(gsub("SignificantGenePairs.csv","",FileNameGenePairs),"OverlapCutOff",toupper(as.character(SignificanceLevel)),"_Edgelist_Biclusters.csv",sep = ""))
  
  write.table(Edgelist.df,file = File.Name,row.names = F,col.names = T,sep = ",")
  
}
