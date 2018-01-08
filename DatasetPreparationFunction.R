#Function to clean input file - takes in gene expression input file with genes along rows and samples along columns with no NAs
DatasetPreparationFunction <- function(InputFileName)
{
  library(data.table)
  df <- fread(InputFileName)
  Gene.IDs <- df[,1]
  Exprs.Matrix <- as.matrix(df[,-1])
  Sample.IDs <- colnames(df)[-1]
  Sums.of.Rows <- apply(Exprs.Matrix,1,sum)
  if (length(which(Sums.of.Rows == 0)) != 0)
  {
    Exprs.Matrix <- Exprs.Matrix[-which(Sums.of.Rows == 0),]
    Non.Zero.Genes <- Gene.IDs[-which(Sums.of.Rows == 0)]
  }
  else
  {
    Non.Zero.Genes <- Gene.IDs
  }
  
  Non.Zero.Genes.df <- data.frame(Non.Zero.Genes)
  colnames(Non.Zero.Genes.df) <- c("Gene.ID")
  Name.EndStr <- substr(InputFileName,nchar(InputFileName)-3,nchar(InputFileName))
  File.Name <- as.character(paste(gsub(Name.EndStr,"",InputFileName),"_GeneNames.csv",sep = ""))
  write.table(Non.Zero.Genes.df,file = "GeneNames.csv",row.names = F,col.names = T)
  Cleaned.df <- data.frame(Non.Zero.Genes,Exprs.Matrix)
  colnames(Cleaned.df) <- c("Gene.ID",Sample.IDs)
  File.Name <- as.character(paste(gsub(Name.EndStr,"",InputFileName),"_Cleaned.csv",sep = ""))
  write.table(Cleaned.df,file = File.Name,row.names = F,col.names = T,sep = ",")
}