
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' @import utils
#' @import stats

#' @title  DataCleaning
#' @description  This function performs basic cleaning of gene expression data sets - removes genes with zero counts across all samples as well as genes with NAs in some samples.
#' @export
#' @param File A character variable. Specifies the name of the .csv or .txt file that contains the gene expression data. Make sure it is in the genes along rows and samples along columns format.
#' @return A data frame containing the gene expression data after basic cleaning. In addition, it generates two .csv files in the working directory - one contains the names of the genes in the cleaned gene expression data set, and the other contains the entire gene expression data after cleaning.
#' @examples
#' \dontrun{
#' DataCleaning(File = "RPGenes.csv")
#' }

DataCleaning <- function(File)
{
  df <- data.table::fread(File)
  Gene.IDs <- df[,1]
  Exprs.Matrix <- as.matrix(df[,-1])
  Sample.IDs <- colnames(df)[-1]
  Sums.of.Rows <- rowSums(Exprs.Matrix)
  if (length(which(Sums.of.Rows == 0 | is.na(Sums.of.Rows) == T)) != 0) {
    Exprs.Matrix <- Exprs.Matrix[-which(Sums.of.Rows == 0 | is.na(Sums.of.Rows) == T),]
    Non.Zero.Genes <- Gene.IDs[-which(Sums.of.Rows == 0 | is.na(Sums.of.Rows) == T)]
  }
  else {
    Non.Zero.Genes <- Gene.IDs
  }

  Non.Zero.Genes.df <- data.frame(Non.Zero.Genes)
  colnames(Non.Zero.Genes.df) <- c("Gene.ID")
  Name.EndStr <- substr(File,nchar(File)-3,nchar(File))
  File.Name <- as.character(paste(gsub(Name.EndStr,"",File),"_GeneNames.csv",sep = ""))
  write.table(Non.Zero.Genes.df,file = File.Name,row.names = F,col.names = T)
  Cleaned.df <- data.frame(Non.Zero.Genes,Exprs.Matrix)
  colnames(Cleaned.df) <- c("Gene.ID",Sample.IDs)
  File.Name <- as.character(paste(gsub(Name.EndStr,"",File),"_Cleaned.csv",sep = ""))
  write.table(Cleaned.df,file = File.Name,row.names = F,col.names = T,sep = ",")
  return(Cleaned.df)
}


#' @title  GenePairs
#' @description  This function finds gene-pairs that share more than a specified proportion of samples between their percentile sets.
#' @export
#' @param File A character variable. Specifies the name of the .csv or .txt file that contains the normalized counts.
#' @param PercSetSize A numeric variable. Specifies the percentage of samples that should be in the percentile sets (strictly greater than 0 and less than 40).
#' @param JcdInd A numeric variable. Specifies the minimum Jaccard Index for the overlap between the percentile sets of a given gene-pair.
#' @param highORlow A character variable. Specifies whether the percentile sets correspond to the highest expression samples ("h") or the lowest expression samples ("l").
#' @param SampleFilter A logical variable. If TRUE, filters out samples over-represented in percentile sets. Default is FALSE.
#' @return A data frame containing the gene-pairs whose Jaccard indices are greater than the specified threshold (JcdInd). Instead of gene symbols their serial numbers in the input gene expression file are used in order to save space. In addition, this function generates 2 .csv files in the working directory - one contains the gene-pairs and their Jaccard indices, while the other contains the binary matrix (genes along rows, samples along columns) in which the presence of a sample in the percentile set of each gene is indicated by a 1. These 2 files are needed as inputs for the \link[TuBA]{Biclustering} function.
#' @examples
#' \dontrun{
#' # For high expression
#' GenePairs(File = "RPGenes.csv",PercSetSize = 5,JcdInd = 0.2,highORlow = "h")
#' # For low expression
#' GenePairs(File = "RPGenes.csv",PercSetSize = 5,JcdInd = 0.2,highORlow = "l")
#' }
GenePairs <- function(File,PercSetSize,JcdInd,highORlow,SampleFilter = NULL)
{
  df <- data.table::fread(File)
  colnames(df) <- c("Gene.ID",colnames(df)[-1])
  Gene.Names <- as.character(df$Gene.ID)
  Relevant.Exprs.Matrix <- as.matrix(df[,-1])
  Sample.IDs <- colnames(df)[-1]
  if (round(PercSetSize) <= 0 | PercSetSize > 40){
    stop("PercSetSize has to be a numeric value strictly greater than 0 and less than 40.")
  }

  CutOffPerc <- round(PercSetSize)/100

  if (JcdInd <= 0 | JcdInd > 1){
    stop("JcdInd has to be a numeric value strictly between 0 and 1.")
  }

  if (is.null(SampleFilter)){
    SampleFilter = F
  }

  if (is.logical(SampleFilter) == F){
    stop("SampleFilter is a logical variable and can only take T or F as valid inputs.")
  }

  if (highORlow == "h" | highORlow == "H" | highORlow == "high" | highORlow == "High" | highORlow == "Hi" | highORlow == "HI" | highORlow == "hi" | highORlow == "HIGH"){

    Threshold <- 1 - CutOffPerc

    ZeroProp.Per.Feature <- apply(Relevant.Exprs.Matrix,1,function(x) mean(x==0))

    Relevant.Features <- which(ZeroProp.Per.Feature < Threshold)

    if (length(Relevant.Features) != 0){
      Gene.Names <- Gene.Names[Relevant.Features]
      Relevant.Exprs.Matrix <- Relevant.Exprs.Matrix[Relevant.Features,]
    }

    rm(df)

    #Find p-value for Jaccard index
    No.of.Samples.In.Percentile.Set <- ceiling(ncol(Relevant.Exprs.Matrix)*CutOffPerc)

    No.of.Samples.In.Intersect <- ceiling(2*JcdInd*No.of.Samples.In.Percentile.Set/(1+JcdInd))

    p.val.Jaccard <- sum(dhyper(No.of.Samples.In.Intersect:No.of.Samples.In.Percentile.Set,No.of.Samples.In.Percentile.Set,ncol(Relevant.Exprs.Matrix)-No.of.Samples.In.Percentile.Set,No.of.Samples.In.Percentile.Set))

    #Matrix that contains the order of the samples in increasing order of expression for all genes
    Sample.Order.Increasing.Matrix <- t(apply(Relevant.Exprs.Matrix,1,order))

    Start.Index <- ncol(Relevant.Exprs.Matrix) - (ceiling(CutOffPerc*ncol(Relevant.Exprs.Matrix))) +1

    #List of outlier samples for every gene
    List.OutlierSamples.SerNos <- vector(mode="list")
    for (i in 1:nrow(Relevant.Exprs.Matrix))
    {
      List.OutlierSamples.SerNos[[i]] <- Sample.Order.Increasing.Matrix[i,Start.Index:ncol(Relevant.Exprs.Matrix)]
    }

    #Binary matrix with 1 for samples that are in the percentile set for any given gene
    Binary.Matrix.For.Genes.Outliers <- matrix(0,nrow = nrow(Sample.Order.Increasing.Matrix),ncol = length(Sample.IDs))

    for (i in 1:nrow(Binary.Matrix.For.Genes.Outliers))
    {
      Binary.Matrix.For.Genes.Outliers[i,List.OutlierSamples.SerNos[[i]]] <- 1
    }

    #Find the frequencies for the samples
    Sample.Frequencies <- colSums(Binary.Matrix.For.Genes.Outliers)

    #Identify samples that only show up in one percentile set
    NonInformativeSamples <- which(Sample.Frequencies == 1)
    if (length(NonInformativeSamples) != 0){
      Binary.Matrix.For.Genes.Outliers <- Binary.Matrix.For.Genes.Outliers[,-NonInformativeSamples]
      Sample.IDs <- Sample.IDs[-NonInformativeSamples]
    }

    if (SampleFilter == T){
      #Find the threshold for maximum frequency based on percentile set size
      MaxThreshold <- max(Sample.Frequencies)
      n <- choose(n = nrow(Binary.Matrix.For.Genes.Outliers),k = 2)
      a <- choose(n = MaxThreshold,k = 2)
      p.val <- a/n
      while (p.val > CutOffPerc) {
        MaxThreshold <- MaxThreshold -1
        a <- choose(n = MaxThreshold,k = 2)
        p.val <- a/n
      }

      NonInformativeSamples <- which(Sample.Frequencies > MaxThreshold)
      if (length(NonInformativeSamples) != 0){
        Binary.Matrix.For.Genes.Outliers <- Binary.Matrix.For.Genes.Outliers[,-NonInformativeSamples]
        Sample.IDs <- Sample.IDs[-NonInformativeSamples]
      }
    }


    Genes.Samples.Binary.df <- data.frame(Gene.Names,Binary.Matrix.For.Genes.Outliers)
    colnames(Genes.Samples.Binary.df) <- c("Gene.ID",Sample.IDs)

    File.Name <- as.character(paste(gsub(".csv","",File),"_H",as.character(CutOffPerc),"_JcdInd",JcdInd,"_GenesSamples_BinaryMatrix.csv",sep = ""))

    write.table(Genes.Samples.Binary.df,file = File.Name,row.names = F,col.names = T,sep = ",")

    rm(Genes.Samples.Binary.df)

    #Outlier overlaps matrix
    SampleOverlaps.Matrix <- matrix(0,nrow = nrow(Relevant.Exprs.Matrix),ncol = nrow(Relevant.Exprs.Matrix))
    SampleOverlaps.Matrix <- Binary.Matrix.For.Genes.Outliers %*% t(Binary.Matrix.For.Genes.Outliers)

    rm(Sample.Order.Increasing.Matrix)

    SampleUnion.Matrix <- matrix(ceiling(2*CutOffPerc*ncol(Relevant.Exprs.Matrix)),nrow = nrow(SampleOverlaps.Matrix),ncol = ncol(SampleOverlaps.Matrix))
    SampleUnion.Matrix <- SampleUnion.Matrix - SampleOverlaps.Matrix

    Jaccard.Dist.Mat <- SampleOverlaps.Matrix/SampleUnion.Matrix
    Jaccard.Dist.Mat[lower.tri(Jaccard.Dist.Mat,diag = T)] <- -1

    #Find relevant gene-pairs

    Relevant.Gene.Pairs <- which(Jaccard.Dist.Mat > JcdInd,arr.ind = T)

    Gene.Pairs.Col1 <- Relevant.Gene.Pairs[,1]

    Gene.Pairs.Col2 <- Relevant.Gene.Pairs[,2]

    Gene.Pairs.df <- data.frame(Gene.Pairs.Col1,Gene.Pairs.Col2,Jaccard.Dist.Mat[Relevant.Gene.Pairs])

    colnames(Gene.Pairs.df) <- c("Gene.1","Gene.2","Jaccard.Index")

    File.Name <- as.character(paste(gsub(".csv","",File),"_H",as.character(CutOffPerc),"_JcdInd",JcdInd,"_GenePairs.csv",sep = ""))

    write.table(Gene.Pairs.df,file = File.Name, row.names = F,col.names = T,sep = ",")

    if (length(Gene.Pairs.df$Gene.1) != 0)
      message(paste0("Successfully generated .csv file containing ",length(Gene.Pairs.df$Gene.1)," gene-pairs (edges)."))


  } else if (highORlow == "l" | highORlow == "L" | highORlow == "low" | highORlow == "Low" | highORlow == "Lo" | highORlow == "lo" | highORlow == "LO" | highORlow == "LOW") {

    Threshold <- CutOffPerc

    ZeroProp.Per.Feature <- apply(Relevant.Exprs.Matrix,1,function(x) mean(x==0))

    Relevant.Features <- which(ZeroProp.Per.Feature <= Threshold)

    if (length(Relevant.Features) != 0){
      Gene.Names <- Gene.Names[Relevant.Features]
      Relevant.Exprs.Matrix <- Relevant.Exprs.Matrix[Relevant.Features,]
    }

    rm(df)

    #Matrix that contains the order of the samples in increasing order of expression for all genes
    Sample.Order.Increasing.Matrix <- t(apply(Relevant.Exprs.Matrix,1,order))

    Start.Index <- 1

    #List of outlier samples for every gene
    List.OutlierSamples.SerNos <- vector(mode="list")
    for (i in 1:nrow(Relevant.Exprs.Matrix))
    {
      List.OutlierSamples.SerNos[[i]] <- Sample.Order.Increasing.Matrix[i,Start.Index:ceiling(CutOffPerc*ncol(Relevant.Exprs.Matrix))]
    }

    #Binary matrix with 1 for samples that are in the percentile set for any given gene
    Binary.Matrix.For.Genes.Outliers <- matrix(0,nrow = nrow(Sample.Order.Increasing.Matrix),ncol = length(Sample.IDs))

    for (i in 1:nrow(Binary.Matrix.For.Genes.Outliers))
    {
      Binary.Matrix.For.Genes.Outliers[i,List.OutlierSamples.SerNos[[i]]] <- 1
    }

    #Find the frequencies for the samples
    Sample.Frequencies <- colSums(Binary.Matrix.For.Genes.Outliers)

    #Identify samples that only show up in one percentile set
    NonInformativeSamples <- which(Sample.Frequencies == 1)
    if (length(NonInformativeSamples) != 0){
      Binary.Matrix.For.Genes.Outliers <- Binary.Matrix.For.Genes.Outliers[,-NonInformativeSamples]
      Sample.IDs <- Sample.IDs[-NonInformativeSamples]
    }

    if (SampleFilter == T){
      #Find the threshold for maximum frequency based on percentile set size
      MaxThreshold <- max(Sample.Frequencies)
      n <- choose(n = nrow(Binary.Matrix.For.Genes.Outliers),k = 2)
      a <- choose(n = MaxThreshold,k = 2)
      p.val <- a/n
      while (p.val > CutOffPerc) {
        MaxThreshold <- MaxThreshold -1
        a <- choose(n = MaxThreshold,k = 2)
        p.val <- a/n
      }

      NonInformativeSamples <- which(Sample.Frequencies > MaxThreshold)
      if (length(NonInformativeSamples) != 0){
        Binary.Matrix.For.Genes.Outliers <- Binary.Matrix.For.Genes.Outliers[,-NonInformativeSamples]
        Sample.IDs <- Sample.IDs[-NonInformativeSamples]
      }
    }

    Genes.Samples.Binary.df <- data.frame(Gene.Names,Binary.Matrix.For.Genes.Outliers)
    colnames(Genes.Samples.Binary.df) <- c("Gene.ID",Sample.IDs)

    File.Name <- as.character(paste(gsub(".csv","",File),"_L",as.character(CutOffPerc),"_JcdInd",JcdInd,"_GenesSamples_BinaryMatrix.csv",sep = ""))

    write.table(Genes.Samples.Binary.df,file = File.Name,row.names = F,col.names = T,sep = ",")

    rm(Genes.Samples.Binary.df)

    #Outlier overlaps matrix
    SampleOverlaps.Matrix <- matrix(0,nrow = nrow(Relevant.Exprs.Matrix),ncol = nrow(Relevant.Exprs.Matrix))
    SampleOverlaps.Matrix <- Binary.Matrix.For.Genes.Outliers %*% t(Binary.Matrix.For.Genes.Outliers)

    rm(Sample.Order.Increasing.Matrix)

    SampleUnion.Matrix <- matrix(ceiling(2*CutOffPerc*ncol(Relevant.Exprs.Matrix)),nrow = nrow(SampleOverlaps.Matrix),ncol = ncol(SampleOverlaps.Matrix))
    SampleUnion.Matrix <- SampleUnion.Matrix - SampleOverlaps.Matrix

    Jaccard.Dist.Mat <- SampleOverlaps.Matrix/SampleUnion.Matrix
    Jaccard.Dist.Mat[lower.tri(Jaccard.Dist.Mat,diag = T)] <- -1

    #Find Jaccard distance between gene-pairs

    Relevant.Gene.Pairs <- which(Jaccard.Dist.Mat > JcdInd,arr.ind = T)

    Gene.Pairs.Col1 <- Relevant.Gene.Pairs[,1]

    Gene.Pairs.Col2 <- Relevant.Gene.Pairs[,2]

    Gene.Pairs.df <- data.frame(Gene.Pairs.Col1,Gene.Pairs.Col2,Jaccard.Dist.Mat[Relevant.Gene.Pairs])

    colnames(Gene.Pairs.df) <- c("Gene.1","Gene.2","Jaccard.Index")

    File.Name <- as.character(paste(gsub(".csv","",File),"_L",as.character(CutOffPerc),"_JcdInd",JcdInd,"_GenePairs.csv",sep = ""))

    write.table(Gene.Pairs.df,file = File.Name, row.names = F,col.names = T,sep = ",")

    if (length(Gene.Pairs.df$Gene.1) != 0)
      message(paste0("Successfully generated .csv file containing ",length(Gene.Pairs.df$Gene.1)," gene-pairs (edges)."))

  }
  return(Gene.Pairs.df)
}

#' @title  Biclustering Function
#' @description  This function finds the biclusters using the graph (gene-pairs file) and the genes-samples binary matrix generated by the \link[TuBA]{GenePairs} function.
#' @export
#' @param VariablePairs A character variable. Specifies the name of the gene-pairs .csv file generated by the \link[TuBA]{GenePairs} function.
#' @param BinaryMatrix A character variable. Specifies the name of the genes-samples binary matrix .csv file generated by the \link[TuBA]{GenePairs} function.
#' @param MinGenes A numeric (integer) variable. Specifies the minimum number of genes that the biclusters must contain.
#' @param MinSamples A numeric (integer) variable. Specifies the minimum number of samples that the biclusters must contain.
#' @param SampleEnrichment A numeric variable - should be greater than 0 and less than or equal to 1. Specifies to what extent the samples should be enriched in the bicluster. Smaller values indicate higher enrichment. Default is 1.
#' @return A data frame containing the sets of genes in the biclusters along with the information about the number of samples in each bicluster. In addition, it generates 3 .csv files - one contains the list of genes in each bicluster along with information about the total number of samples in the bicluster and how many of them are contributed by each gene within the bicluster; the other contains the bicluster-samples binary matrix (biclusters along rows, samples along columns) in which the presence of a sample within a bicluster is indicated by a 1; the third file contains the genes-biclusters-samples matrix (genes along rows, samples along columns) which contains more detailed information about which sample is present for which gene within a given bicluster. The first column in this genes-biclusters-samples file contains the names of genes in the bicluster, and the second column contains the serial numbers of the biclusters which contains these genes; the rest of the file contains the binary matrix.
#' @examples
#' \dontrun{
#' Biclustering(VariablePairs = "RPGenes_H0.05_JcdInd0.2_GenePairs.csv",BinaryMatrix = "RPGenes_H0.05_JcdInd0.2_GenesSamples_BinaryMatrix.csv")
#' }
Biclustering <- function(VariablePairs,BinaryMatrix,MinGenes = NULL,MinSamples = NULL,SampleEnrichment = NULL)
{
  if(is.null(MinGenes)){
    MinGenes <- 3
    message("MinGenes not specified. Minimum of 3 set as default.")
  }

  if(MinGenes < 3){
    MinGenes <- 3
    message("MinGenes cannot be less than 3. Minimum of 3 set as default.")
  }

  if(is.null(MinSamples)){
    MinSamples <- 2
    message("MinSamples not specified. Will seek biclusters containing at least 2 samples.")
  }

  if (is.null(SampleEnrichment)){
    SampleEnrichment <- 1
  }

  #Import data frame that contains the variable pairs found by the significant node pairs function
  Variable.Pairs.df <- data.table::fread(VariablePairs)
  colnames(Variable.Pairs.df) <- c("Variable.1","Variable.2","JcdInd")

  JaccardOverlapProp <- round(min(Variable.Pairs.df$JcdInd),digits = 2)

  #Import data frame that contains the binary matrix of variables and their respective percentile set samples
  Variables.Samples.Binary.df <- data.table::fread(BinaryMatrix)
  colnames(Variables.Samples.Binary.df) <- c("Variable.ID",colnames(Variables.Samples.Binary.df[,-1]))

  #Names of variables
  Variable.Names <- as.character(Variables.Samples.Binary.df$Variable.ID)

  #Binary matrix of variables and percentile set samples
  Matrix.For.Variables.Outliers <- as.matrix(Variables.Samples.Binary.df[,-1])
  Binary.Matrix.For.Variables.Outliers <- matrix(0,nrow = nrow(Matrix.For.Variables.Outliers),ncol = ncol(Matrix.For.Variables.Outliers))
  Binary.Matrix.For.Variables.Outliers[which(Matrix.For.Variables.Outliers == 1 | Matrix.For.Variables.Outliers == 11,arr.ind = T)] <- 1
  colnames(Binary.Matrix.For.Variables.Outliers) <- colnames(Matrix.For.Variables.Outliers)

  #IDs of samples (or conditions)
  Sample.IDs <- colnames(Matrix.For.Variables.Outliers)

  rm(Variables.Samples.Binary.df)

  if(length(Variable.Pairs.df$Variable.1) == 0){
    stop("Please check input file. No variable pairs found.")
  } else {
    message(paste0("There are ",length(Variable.Pairs.df$Variable.1)," edges in the graph."))
  }

  #Column 1 of qualified variable pairs data frame
  Qualified.Variable1 <- Variable.Pairs.df$Variable.1

  #Column 2 of qualified variable pairs data frame
  Qualified.Variable2 <- Variable.Pairs.df$Variable.2

  #Summarize graph info in a table
  Variable.Summary.Info <- table(c(Qualified.Variable1,Qualified.Variable2))

  #All variables in graph
  All.Variables.In.Graph <- as.numeric(names(Variable.Summary.Info))

  #Degrees of variables in graph
  Degrees.Of.Variables <- as.numeric(Variable.Summary.Info)

  rm(Variable.Pairs.df)

  #Prepare adjacency matrix for the graph
  Adjacency.Mat.Variables <- matrix(0,nrow = length(Variable.Names),ncol = length(Variable.Names))

  Adjacency.Mat.Variables[cbind(Qualified.Variable1,Qualified.Variable2)] <- 1

  message("Preparing graph...")

  #Make the adjacency matrix symmetric
  Adjacency.Mat.Variables <- Adjacency.Mat.Variables + t(Adjacency.Mat.Variables)

  #Reduce the adjacency matrix such that it only contains variables that are present in the graph
  Adjacency.Mat.Variables <- Adjacency.Mat.Variables[All.Variables.In.Graph,All.Variables.In.Graph]
  rownames(Adjacency.Mat.Variables) <- All.Variables.In.Graph
  colnames(Adjacency.Mat.Variables) <- All.Variables.In.Graph

  #Vector to convert variable serial numbers to the corresponding serial number in the reduced adjacency matrix
  Variable.Annotation.Conversion.Vec <- vector(mode = "numeric", length = length(Variable.Names))
  Variable.Ser.Nos <- 1:length(Variable.Names)
  Matching.Variables.Indices <- match(All.Variables.In.Graph,Variable.Ser.Nos)
  Variable.Annotation.Conversion.Vec[Matching.Variables.Indices] <- 1:length(All.Variables.In.Graph)

  message("Finding biclusters...")

  if (length(Qualified.Variable1) > 200000)
    message("This may take several minutes due to the large size of the graph.")

  Total.No.of.Edges.In.Unpruned.Graph <- length(Qualified.Variable1)

  #Find the number of triangular cliques associated with each node-pair in graph
  No.of.Nodes.Associated <- vector(mode = "numeric",length = length(Qualified.Variable1))
  Samples.Background.Frequencies <- vector(mode = "numeric",length = length(Sample.IDs))
  for (i in 1:length(Qualified.Variable1))
  {
    TempI1 <- Qualified.Variable1[i]
    TempI2 <- Qualified.Variable2[i]
    Adjacency.Vec1 <- Adjacency.Mat.Variables[Variable.Annotation.Conversion.Vec[TempI1],]
    Adjacency.Vec2 <- Adjacency.Mat.Variables[Variable.Annotation.Conversion.Vec[TempI2],]
    No.of.Nodes.Associated[i] <- length(which(Adjacency.Vec1*Adjacency.Vec2 == 1))
    Sample.Set1 <- Binary.Matrix.For.Variables.Outliers[TempI1,]
    Sample.Set2 <- Binary.Matrix.For.Variables.Outliers[TempI2,]
    Temp.Sample.Ser.Nos <- which(Sample.Set1*Sample.Set2 == 1)
    Samples.Background.Frequencies[Temp.Sample.Ser.Nos] <- Samples.Background.Frequencies[Temp.Sample.Ser.Nos] + 1
  }

  #Non-triangular gene-pairs
  Node.Pairs.For.Filtering <- which(No.of.Nodes.Associated == 0)

  if (length(Node.Pairs.For.Filtering) != 0){
    No.of.Nodes.Associated <- No.of.Nodes.Associated[-Node.Pairs.For.Filtering]

    Qualified.Variable1 <- Qualified.Variable1[-Node.Pairs.For.Filtering]
    Qualified.Variable2 <- Qualified.Variable2[-Node.Pairs.For.Filtering]
  }

  if(max(No.of.Nodes.Associated) == 0){
    stop("No clique of size 3 found in input graph.")
  }

  Decreasing.No.of.Nodes <- order(No.of.Nodes.Associated,decreasing = T)

  #Sort column 1 and column2 based on decreasing order of nodes associated
  Qualified.Variable1 <- Qualified.Variable1[Decreasing.No.of.Nodes]
  Qualified.Variable2 <- Qualified.Variable2[Decreasing.No.of.Nodes]

  #Find potential biclusters
  i <- 1
  Temp.Vec1 <- Qualified.Variable1
  Temp.Vec2 <- Qualified.Variable2
  Nodes.In.Subgraph <- vector(mode = "list")
  while (length(Temp.Vec1) != 0){
    Adjacency.Vec1 <- Adjacency.Mat.Variables[Variable.Annotation.Conversion.Vec[Temp.Vec1[1]],]
    Adjacency.Vec2 <- Adjacency.Mat.Variables[Variable.Annotation.Conversion.Vec[Temp.Vec2[1]],]
    Nodes.In.Subgraph[[i]] <- c(Temp.Vec1[1],Temp.Vec2[1],All.Variables.In.Graph[which(Adjacency.Vec1*Adjacency.Vec2 == 1)])

    #Remove those edges that contain the variables in the dense subgraph
    Edges.To.Be.Removed <- which(Temp.Vec1 %in% Nodes.In.Subgraph[[i]] | Temp.Vec2 %in% Nodes.In.Subgraph[[i]])

    Temp.Vec1 <- Temp.Vec1[-Edges.To.Be.Removed]
    Temp.Vec2 <- Temp.Vec2[-Edges.To.Be.Removed]

    i <- i + 1
  }

  #Reintroduce dense subgraphs back in original graph and add nodes that share edges with at least 2 nodes in dense subgraph
  Nodes.In.Bicluster <- vector(mode = "list",length = length(Nodes.In.Subgraph))
  for (i in 1:length(Nodes.In.Subgraph))
  {
    Col.Ser.Nos.For.Bicluster <- Variable.Annotation.Conversion.Vec[Nodes.In.Subgraph[[i]]]
    Sub.Adj.Mat <- Adjacency.Mat.Variables[,Col.Ser.Nos.For.Bicluster]
    Nodes.In.Bicluster[[i]] <- unique(c(Nodes.In.Subgraph[[i]],as.numeric(rownames(Sub.Adj.Mat)[which(rowSums(Sub.Adj.Mat) >= 2)])))
  }

  Bicluster.Size.Order <- order(unlist(lapply(Nodes.In.Bicluster,length)),decreasing = T)

  Nodes.In.Bicluster <- Nodes.In.Bicluster[Bicluster.Size.Order]

  #Find biclusters that have nodes that are subsets of other biclusters
  Nested.Biclusters <- vector(mode = "numeric")
  for (i in 1:length(Nodes.In.Bicluster))
  {
    TempI <- length(Nodes.In.Bicluster) - (i-1)
    Temp.Bicluster.I <- Nodes.In.Bicluster[[TempI]]
    j <- 1
    Temp.Intersection <- 1
    while(Temp.Intersection != 0 & TempI != j){
      Temp.Bicluster.J <- Nodes.In.Bicluster[[j]]
      Temp.Intersection <- length(Temp.Bicluster.I) - length(intersect(Temp.Bicluster.I,Temp.Bicluster.J))
      if (Temp.Intersection == 0 & TempI != j & length(Temp.Bicluster.J) >= length(Temp.Bicluster.I)){
        Nested.Biclusters <- c(Nested.Biclusters,TempI)
      }
      j <- j + 1
    }
  }

  #Remove biclusters that are nested within larger biclusters
  if (length(Nested.Biclusters) != 0){
    Nodes.In.Bicluster <- Nodes.In.Bicluster[-Nested.Biclusters]
  }

  #Filter out biclusters that have fewer genes than the specified threshold (MinGenes)
  if (length(which(unlist(lapply(Nodes.In.Bicluster,length)) >= MinGenes)) == 0){
    stop("No biclusters with ",MinGenes," genes found. Please choose different parameters.")
  } else {
    Nodes.In.Bicluster <- Nodes.In.Bicluster[which(unlist(lapply(Nodes.In.Bicluster,length)) >= MinGenes)]
  }

  #Find samples preferentially associated with biclusters
  n <- Total.No.of.Edges.In.Unpruned.Graph
  Bicluster.Samples.Matrix <- matrix(0,nrow = length(Nodes.In.Bicluster),ncol = ncol(Binary.Matrix.For.Variables.Outliers))
  Samples.In.Bicluster <- vector(mode = "list",length = length(Nodes.In.Bicluster))
  for (i in 1:length(Nodes.In.Bicluster))
  {
    Temp.Nodes.In.Bicluster <- Nodes.In.Bicluster[[i]]
    Temp.Ser.Nos.Edges.In.Subgraph <- which(Qualified.Variable1 %in% Temp.Nodes.In.Bicluster & Qualified.Variable2 %in% Temp.Nodes.In.Bicluster)

    No.of.Edges.In.Subgraph <- length(Temp.Ser.Nos.Edges.In.Subgraph)

    Bicluster.Samples.Frequencies <- vector(mode = "numeric",length = length(Sample.IDs))
    for (j in 1:length(Temp.Ser.Nos.Edges.In.Subgraph))
    {
      TempJ1 <- Qualified.Variable1[Temp.Ser.Nos.Edges.In.Subgraph[j]]
      TempJ2 <- Qualified.Variable2[Temp.Ser.Nos.Edges.In.Subgraph[j]]
      Temp.Sample.Set1 <- Binary.Matrix.For.Variables.Outliers[TempJ1,]
      Temp.Sample.Set2 <- Binary.Matrix.For.Variables.Outliers[TempJ2,]
      Temp.Sample.Ser.Nos <- which(Temp.Sample.Set1*Temp.Sample.Set2 == 1)
      Bicluster.Samples.Frequencies[Temp.Sample.Ser.Nos] <- Bicluster.Samples.Frequencies[Temp.Sample.Ser.Nos] + 1
    }

    Valid.Sample.Ser.Nos <- which(Bicluster.Samples.Frequencies != 0)
    Bicluster.Samples.Frequencies <- Bicluster.Samples.Frequencies[Valid.Sample.Ser.Nos]
    Order.Bicluster.Samples.Frequencies <- order(Bicluster.Samples.Frequencies,decreasing = T)
    Bicluster.Samples.Frequencies <- Bicluster.Samples.Frequencies[Order.Bicluster.Samples.Frequencies]
    Valid.Sample.Ser.Nos <- Valid.Sample.Ser.Nos[Order.Bicluster.Samples.Frequencies]

    if (SampleEnrichment == 1){
      Sample.Ratios.Per.Bicluster <- vector(mode = "numeric",length = length(Valid.Sample.Ser.Nos))
      for (k in 1:length(Valid.Sample.Ser.Nos))
      {
        t <- Bicluster.Samples.Frequencies[k]
        Sample.Count.In.Graph <- Samples.Background.Frequencies[Valid.Sample.Ser.Nos[k]]

        Sample.Ratios.Per.Bicluster[k] <- (t/Sample.Count.In.Graph)/(No.of.Edges.In.Subgraph/n)
      }

      Temp.ser.nos <- which(Sample.Ratios.Per.Bicluster >= 1)
      Samples.In.Bicluster[[i]] <- Valid.Sample.Ser.Nos[Temp.ser.nos]
      Bicluster.Samples.Matrix[i,Valid.Sample.Ser.Nos[Temp.ser.nos]] <- 1
    } else {
      p.value.sample <- vector(mode = "numeric",length = length(Valid.Sample.Ser.Nos))
      Sample.Ratios.Per.Bicluster <- vector(mode = "numeric",length = length(Valid.Sample.Ser.Nos))
      for (k in 1:length(Valid.Sample.Ser.Nos))
      {
        Temp.Sample.Frequency <- Bicluster.Samples.Frequencies[k]
        t <- Temp.Sample.Frequency

        Sample.Count.In.Graph <- Samples.Background.Frequencies[Valid.Sample.Ser.Nos[k]]
        if (No.of.Edges.In.Subgraph <= Sample.Count.In.Graph){
          b <- No.of.Edges.In.Subgraph
          a <- Sample.Count.In.Graph
        } else {
          a <- No.of.Edges.In.Subgraph
          b <- Sample.Count.In.Graph
        }
        p.value.sample[k] <- sum(dhyper(t:b,a,n-a,b))
        Sample.Ratios.Per.Bicluster[k] <- (t/Sample.Count.In.Graph)/(No.of.Edges.In.Subgraph/n)
      }
      Temp.ser.nos <- which(p.value.sample <= SampleEnrichment & Sample.Ratios.Per.Bicluster >= 1)
      Samples.In.Bicluster[[i]] <- Valid.Sample.Ser.Nos[Temp.ser.nos]
      Bicluster.Samples.Matrix[i,Valid.Sample.Ser.Nos[Temp.ser.nos]] <- 1
    }
  }

  #Find biclusters that contain at least the minimum of samples specified (MinSamples)
  Biclusters.With.Some.Samples <- which(rowSums(Bicluster.Samples.Matrix) >= MinSamples)

  if (length(Biclusters.With.Some.Samples) == 0 & SampleEnrichment != 1){
    stop("No biclusters were found with the given choice of SampleEnrichment. Try with SampleEnrichment = 1 (default)")
  } else if (length(Biclusters.With.Some.Samples) == 0 & SampleEnrichment == 1){
    stop("No bicluster found with at least 3 genes")
  }

  #Filter nodes in biclusters based on the samples found enriched in each bicluster - Approach 1 (Remove nodes based on the overlap of their percentile sets with samples found enriched in the bicluster)
  Nodes.In.Bicluster.Eligible <- Nodes.In.Bicluster[Biclusters.With.Some.Samples]

  Bicluster.Samples.Matrix.Eligible <- Bicluster.Samples.Matrix[Biclusters.With.Some.Samples,]

  Nodes.In.Final.Biclusters <- vector(mode = "list",length = length(Nodes.In.Bicluster.Eligible))
  N.PercentileSet <- sum(Binary.Matrix.For.Variables.Outliers[1,])
  for (i in 1:length(Nodes.In.Bicluster.Eligible))
  {
    Temp.Nodes.In.Bicluster <- Nodes.In.Bicluster.Eligible[[i]]
    JInd <- vector(mode = "numeric",length = length(Temp.Nodes.In.Bicluster))
    for (j in 1:length(Temp.Nodes.In.Bicluster))
    {
      Samples.In.Percentile.SetJ <- which(Binary.Matrix.For.Variables.Outliers[Temp.Nodes.In.Bicluster[j],] == 1)
      Intersecting.Samples <- intersect(which(Bicluster.Samples.Matrix[i,] == 1),Samples.In.Percentile.SetJ)
      JInd[j] <- length(Intersecting.Samples)/(2*N.PercentileSet -  length(Intersecting.Samples))
    }
    Temp.Filter.Index <- which(JInd < JaccardOverlapProp)
    if (length(Temp.Filter.Index != 0))
    {
      Nodes.In.Final.Biclusters[[i]] <- Nodes.In.Bicluster.Eligible[[i]][-Temp.Filter.Index]
    } else {
      Nodes.In.Final.Biclusters[[i]] <- Nodes.In.Bicluster.Eligible[[i]]
    }
  }

  #Remove biclusters that have fewer nodes than MinGenes
  Empty.Final.Biclusters <- which(unlist(lapply(Nodes.In.Final.Biclusters,length)) <= MinGenes-1)

  if (length(Empty.Final.Biclusters) != 0){
    Nodes.In.Final.Biclusters <- Nodes.In.Final.Biclusters[-Empty.Final.Biclusters]
    Bicluster.Samples.Matrix.Eligible <- Bicluster.Samples.Matrix.Eligible[-Empty.Final.Biclusters,]
  }

  Nested.Biclusters <- vector(mode = "numeric")
  if (length(Nodes.In.Final.Biclusters) != 0){
    for (i in 1:length(Nodes.In.Final.Biclusters))
    {
      TempI <- length(Nodes.In.Final.Biclusters) - (i-1)
      Temp.Bicluster.I <- Nodes.In.Final.Biclusters[[TempI]]
      j <- 1
      Temp.Intersection <- 1
      while(Temp.Intersection != 0 & TempI != j){
        Temp.Bicluster.J <- Nodes.In.Final.Biclusters[[j]]
        Temp.Intersection <- length(Temp.Bicluster.I) - length(intersect(Temp.Bicluster.I,Temp.Bicluster.J))
        if (Temp.Intersection == 0 & TempI != j & length(Temp.Bicluster.J) >= length(Temp.Bicluster.I)){
          Nested.Biclusters <- c(Nested.Biclusters,TempI)
        }
        j <- j + 1
      }
    }
  }

  if (length(Nested.Biclusters) != 0){
    Nodes.In.Final.Biclusters <- Nodes.In.Final.Biclusters[-Nested.Biclusters]
    Bicluster.Samples.Matrix.Eligible <- Bicluster.Samples.Matrix.Eligible[-Nested.Biclusters,]
  }

  #Rearrange biclusters based on the number of nodes they contain (largest to smallest)
  Bicluster.Size.Order <- order(unlist(lapply(Nodes.In.Final.Biclusters,length)),decreasing = T)

  if (length(Bicluster.Size.Order) > 1){
    Nodes.In.Final.Biclusters <- Nodes.In.Final.Biclusters[Bicluster.Size.Order]
    Bicluster.Samples.Matrix.Eligible <- Bicluster.Samples.Matrix.Eligible[Bicluster.Size.Order,]
  } else {
    Bicluster.Samples.Matrix.Eligible <- t(as.matrix(Bicluster.Samples.Matrix.Eligible))
  }

  if (length(dim(Bicluster.Samples.Matrix.Eligible)) == 0){
    MinBiclusterSamples <- sum(Bicluster.Samples.Matrix.Eligible)
    No.of.Samples <- length(Bicluster.Samples.Matrix.Eligible)
  } else if (length(rowSums(Bicluster.Samples.Matrix.Eligible)) != 0){
    MinBiclusterSamples <- min(rowSums(Bicluster.Samples.Matrix.Eligible))
    No.of.Samples <- ncol(Bicluster.Samples.Matrix.Eligible)
  } else {
    MinBiclusterSamples <- 0
  }

  Temp.No.of.Nodes.In.Bicluster <- unlist(lapply(Nodes.In.Final.Biclusters,length))
  if (length(Temp.No.of.Nodes.In.Bicluster) != 0){
    MinBiclusterGenes <- min(Temp.No.of.Nodes.In.Bicluster)
  } else {
    MinBiclusterGenes <- 0
  }

  OriginalMinSamples <- MinSamples
  if (MinBiclusterSamples > MinSamples)
    MinSamples <- MinBiclusterSamples

  if (OriginalMinSamples != 2){
    message(paste0("Found ",length(Nodes.In.Final.Biclusters)," biclusters with at least ",MinGenes," genes and ",OriginalMinSamples," samples."))
  } else {
    message(paste0("Found ",length(Nodes.In.Final.Biclusters)," biclusters with at least ",MinGenes," genes and ",MinSamples," samples."))
  }

  #Prepare the output files
  if (length(Nodes.In.Final.Biclusters) != 0){

    Nodes.Biclusters.Info.df <- data.frame(Variable.Names[unlist(Nodes.In.Final.Biclusters)],rep(1:length(Nodes.In.Final.Biclusters),unlist(lapply(Nodes.In.Final.Biclusters,length))))
    colnames(Nodes.Biclusters.Info.df) <- c("Gene.ID","Bicluster.No")

    No.of.Samples.In.Bicluster <- vector(mode = "numeric",length = length(Nodes.Biclusters.Info.df$Gene.ID))
    Samples.Per.Gene.In.Bicluster <- vector(mode = "list",length = length(Nodes.Biclusters.Info.df$Gene.ID))
    No.of.Samples.Per.Gene.In.Bicluster <- vector(mode = "numeric",length = length(Nodes.Biclusters.Info.df$Gene.ID))
    Proportion.of.Samples.In.Bicluster <- vector(mode = "numeric",length = length(Nodes.Biclusters.Info.df$Gene.ID))
    Genes.Bicluster.Samples.Matrix <- matrix(0,nrow = length(Nodes.Biclusters.Info.df$Gene.ID),ncol = No.of.Samples)
    for (i in 1:length(Nodes.Biclusters.Info.df$Gene.ID))
    {
      Gene.Ser.No <- which(Variable.Names == Nodes.Biclusters.Info.df$Gene.ID[i])
      Bicluster.No <- Nodes.Biclusters.Info.df$Bicluster.No[i]
      if (length(dim(Bicluster.Samples.Matrix.Eligible)) == 0){
        No.of.Samples.In.Bicluster[i] <- length(which(Bicluster.Samples.Matrix.Eligible == 1))
        Samples.Per.Gene.In.Bicluster[[i]] <- intersect(which(Binary.Matrix.For.Variables.Outliers[Gene.Ser.No,] == 1),which(Bicluster.Samples.Matrix.Eligible == 1))
      } else {
        No.of.Samples.In.Bicluster[i] <- length(which(Bicluster.Samples.Matrix.Eligible[Bicluster.No,] == 1))
        Samples.Per.Gene.In.Bicluster[[i]] <- intersect(which(Binary.Matrix.For.Variables.Outliers[Gene.Ser.No,] == 1),which(Bicluster.Samples.Matrix.Eligible[Bicluster.No,] == 1))
      }

      No.of.Samples.Per.Gene.In.Bicluster[i] <- length(Samples.Per.Gene.In.Bicluster[[i]])
      Proportion.of.Samples.In.Bicluster[i] <- No.of.Samples.Per.Gene.In.Bicluster[i]/No.of.Samples.In.Bicluster[i]
      Genes.Bicluster.Samples.Matrix[i,Samples.Per.Gene.In.Bicluster[[i]]] <- 1
    }

    Nodes.Biclusters.Info.df$Samples.In.Bicluster <- No.of.Samples.In.Bicluster
    Nodes.Biclusters.Info.df$Samples.Per.Gene.In.Bicluster <- No.of.Samples.Per.Gene.In.Bicluster
    Nodes.Biclusters.Info.df$Proportion.of.Samples <- Proportion.of.Samples.In.Bicluster

    File.Name <- paste0(substr(VariablePairs,1,nchar(VariablePairs)-24),"_MinGenes",MinBiclusterGenes,"_MinSamples",MinBiclusterSamples,"_GenesInBiclusters.csv")

    write.table(Nodes.Biclusters.Info.df,file = File.Name,row.names = F,col.names = T,sep = ",")

    message("Successfully generated .csv file containing the list of genes in biclusters.")

    Bicluster.Samples.df <- data.frame(paste0("Bicluster.",1:length(Nodes.In.Final.Biclusters)),Bicluster.Samples.Matrix.Eligible)
    colnames(Bicluster.Samples.df) <- c("Bicluster.No",colnames(Binary.Matrix.For.Variables.Outliers))

    File.Name <- paste0(substr(VariablePairs,1,nchar(VariablePairs)-24),"_MinGenes",MinBiclusterGenes,"_MinSamples",MinBiclusterSamples,"_BiclusterSamplesMatrix.csv")

    write.table(Bicluster.Samples.df,file = File.Name,row.names = F,col.names = T,sep = ",")

    if (nrow(Bicluster.Samples.Matrix.Eligible) != 0)
      message("Successfully generated .csv file containing the bicluster-samples binary matrix.")

    Genes.Bicluster.Samples.df <- data.frame(Nodes.Biclusters.Info.df$Gene.ID,Nodes.Biclusters.Info.df$Bicluster.No,Genes.Bicluster.Samples.Matrix)
    colnames(Genes.Bicluster.Samples.df) <- c("Gene.ID","Bicluster.No",colnames(Binary.Matrix.For.Variables.Outliers))

    File.Name <- paste0(substr(VariablePairs,1,nchar(VariablePairs)-24),"_MinGenes",MinBiclusterGenes,"_MinSamples",MinBiclusterSamples,"_GenesBiclusterSamplesMatrix.csv")

    write.table(Genes.Bicluster.Samples.df,file = File.Name,row.names = F,col.names = T,sep = ",")

    if (nrow(Bicluster.Samples.Matrix.Eligible) != 0)
      message("Successfully generated .csv file containing the genes-bicluster-samples binary matrix.")

  }
  return(Nodes.Biclusters.Info.df)
}



