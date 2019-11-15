GenesSamplesBiclusterAssociations <- function(InputMatrixFileName,BiclusterSamplesMatrixFileName,GenesBiclustersInfoFileName,highORlow){

  library(data.table)
  if (highORlow == "high" | highORlow == "HIGH" | highORlow == "h" | highORlow == "H" | highORlow == "Hi" | highORlow == "HI"){
    #Expression Matrix

    Exp.df <- fread(InputMatrixFileName)

    Exp.Matrix <- as.matrix(Exp.df[,-1])

    Gene.IDs.Vec <- Exp.df$Gene.ID

    #Matrix that contains the order of the samples in increasing order of expression for all genes

    Sample.Order.Increasing.Matrix <- t(apply(Exp.Matrix,1,order))

    #Read in bicluster sample file

    Bicluster.Sample.df <- fread(BiclusterSamplesMatrixFileName)

    #Matrix with sample membership in biclusters

    Bicluster.Sample.Matrix <- as.matrix(Bicluster.Sample.df[,-1])

    Sample.IDs <- colnames(Bicluster.Sample.Matrix)

    #Read in gene degrees bicluster file

    Gene.Degrees.Bicluster.df <- fread(GenesBiclustersInfoFileName)

    #Number of samples in dataset

    No.of.Samples <- ncol(Bicluster.Sample.Matrix)

    #Serial numbers of samples

    Ser.Nos.Of.Samples <- seq(1,No.of.Samples,1)

    #Ser.Nos of Non-Zero Samples per gene
    List.of.Ser.Nos.of.Non.Zero.Samples <- vector(mode = "list",length = length(Gene.IDs.Vec))
    for (i in 1:length(Gene.IDs.Vec))
    {
      List.of.Ser.Nos.of.Non.Zero.Samples[[i]] <- as.numeric(which(Exp.Matrix[i,] !=0))
    }

    Lengths.Non.Zero.Samples <- unlist(lapply(List.of.Ser.Nos.of.Non.Zero.Samples,length))


    #Code to find enrichments in top X samples for each gene in biclusters

    if (length(which(Lengths.Non.Zero.Samples < No.of.Samples)) != 0){
      List.Of.Enrichments.Per.Gene.Per.Bicluster <- vector(mode = "list",length = nrow(Bicluster.Sample.Matrix))
      for (i in 1:nrow(Bicluster.Sample.Matrix))
      {
        Ser.Nos.Samples.In.Bicluster <- which(Bicluster.Sample.Matrix[i,] == 1)
    
        Sample.IDs.In.Bicluster <- Sample.IDs[Ser.Nos.Samples.In.Bicluster]
    
        Start.Index <- ncol(Exp.Matrix) - length(Ser.Nos.Samples.In.Bicluster) +1
    
        Gene.IDs.In.Bicluster <- Gene.Degrees.Bicluster.df$Gene.ID[which(Gene.Degrees.Bicluster.df$Bicluster.No == i)]
        List.Containing.Sample.Order <- vector(mode = "list",length = length(Gene.IDs.In.Bicluster))
    
        p.val.Overlaps <- vector(mode = "numeric",length = length(Gene.IDs.In.Bicluster))
        for (j in 1:length(Gene.IDs.In.Bicluster))
        {
          if (length(which(Gene.IDs.Vec == Gene.IDs.In.Bicluster[j])) > 1){
            TempJ <- which(Gene.IDs.Vec == Gene.IDs.In.Bicluster[j])[2]
          } else {
            TempJ <- which(Gene.IDs.Vec == Gene.IDs.In.Bicluster[j])
          }
          Sample.IDs.For.GeneJ <- Sample.IDs[Sample.Order.Increasing.Matrix[TempJ,Start.Index:ncol(Exp.Matrix)]]
      
          t <- length(intersect(Sample.IDs.For.GeneJ,Sample.IDs.In.Bicluster))
          n <- length(Sample.IDs)
          a <- length(Sample.IDs.In.Bicluster)
      
          if (Lengths.Non.Zero.Samples[TempJ] < length(Sample.IDs.In.Bicluster)){
            b <- Lengths.Non.Zero.Samples[TempJ]
          } else {
            b <- length(Sample.IDs.For.GeneJ)
          }
      
          p.val.Overlaps[j] <- sum(dhyper(t:b,a,n-a,b))
        }
        List.Of.Enrichments.Per.Gene.Per.Bicluster[[i]] <- p.val.Overlaps 
      }
    } else {
      List.Of.Enrichments.Per.Gene.Per.Bicluster <- vector(mode = "list",length = nrow(Bicluster.Sample.Matrix))
      for (i in 1:nrow(Bicluster.Sample.Matrix))
      {
        Ser.Nos.Samples.In.Bicluster <- which(Bicluster.Sample.Matrix[i,] == 1)
    
        Sample.IDs.In.Bicluster <- Sample.IDs[Ser.Nos.Samples.In.Bicluster]
    
        Start.Index <- ncol(Exp.Matrix) - length(Ser.Nos.Samples.In.Bicluster) +1
    
        Gene.IDs.In.Bicluster <- Gene.Degrees.Bicluster.df$Gene.ID[which(Gene.Degrees.Bicluster.df$Bicluster.No == i)]
        List.Containing.Sample.Order <- vector(mode = "list",length = length(Gene.IDs.In.Bicluster))
    
        p.val.Overlaps <- vector(mode = "numeric",length = length(Gene.IDs.In.Bicluster))
        for (j in 1:length(Gene.IDs.In.Bicluster))
        {
          TempJ <- which(Gene.IDs.Vec == Gene.IDs.In.Bicluster[j])
          Sample.IDs.For.GeneJ <- Sample.IDs[Sample.Order.Increasing.Matrix[TempJ,Start.Index:ncol(Exp.Matrix)]]
      
          t <- length(intersect(Sample.IDs.For.GeneJ,Sample.IDs.In.Bicluster))
          n <- length(Sample.IDs)
          a <- length(Sample.IDs.In.Bicluster)
          b <- length(Sample.IDs.For.GeneJ)
      
          p.val.Overlaps[j] <- sum(dhyper(t:b,a,n-a,b))
        }
        List.Of.Enrichments.Per.Gene.Per.Bicluster[[i]] <- p.val.Overlaps 
      }
  
    }

    Vec.With.Enrichments.From.List <- unlist(List.Of.Enrichments.Per.Gene.Per.Bicluster)

    FDR.Vec.With.Enrichments.From.List <- p.adjust(Vec.With.Enrichments.From.List,method = "BH")

    New.Gene.Degrees.Bicluster.df <- data.frame(Gene.Degrees.Bicluster.df,Vec.With.Enrichments.From.List,FDR.Vec.With.Enrichments.From.List)
    colnames(New.Gene.Degrees.Bicluster.df) <- c(colnames(Gene.Degrees.Bicluster.df),"p.val.Enrichment.of.Top.Samples.In.Bicluster","FDR.Enrichment.of.Top.Samples.In.Bicluster")

    FileName <- paste0(gsub(".csv","",GenesBiclustersInfoFileName),"_WithAssociationPValue.csv")
    
    write.table(New.Gene.Degrees.Bicluster.df,file = FileName,row.names = F,col.names = T,sep = ",")
    
    if (length(which(Lengths.Non.Zero.Samples < No.of.Samples)) != 0){
      List.With.Genes.Membership.Per.Sample.Per.Bicluster <- vector(mode = "list",length = nrow(Bicluster.Sample.Matrix))
      for (i in 1:nrow(Bicluster.Sample.Matrix))
      {
        Ser.Nos.Samples.In.Bicluster <- which(Bicluster.Sample.Matrix[i,] == 1)
        
        Sample.IDs.In.Bicluster <- Sample.IDs[Ser.Nos.Samples.In.Bicluster]
        
        Start.Index <- ncol(Exp.Matrix) - length(Ser.Nos.Samples.In.Bicluster) +1
        
        Gene.IDs.In.Bicluster <- Gene.Degrees.Bicluster.df$Gene.ID[which(Gene.Degrees.Bicluster.df$Bicluster.No == i)]
        
        
        List.Genes.With.SampleK <- vector(mode = "list",length = length(Sample.IDs.In.Bicluster))
        for (k in 1:length(Sample.IDs.In.Bicluster))
        {
          Genes.With.SampleK <- vector(mode = "character")
          for (j in 1:length(Gene.IDs.In.Bicluster))
          {
            if (length(which(Gene.IDs.Vec == Gene.IDs.In.Bicluster[j])) > 1){
              TempJ <- which(Gene.IDs.Vec == Gene.IDs.In.Bicluster[j])[2]
            } else {
              TempJ <- which(Gene.IDs.Vec == Gene.IDs.In.Bicluster[j])
            }
            
            if (Lengths.Non.Zero.Samples[TempJ] < length(Sample.IDs.In.Bicluster)){
              Sample.IDs.For.GeneJ <- Sample.IDs[which(Exp.Matrix[TempJ,] != 0)]
            } else {
              Sample.IDs.For.GeneJ <- Sample.IDs[Sample.Order.Increasing.Matrix[TempJ,Start.Index:ncol(Exp.Matrix)]]
            }
            
            Sample.ID.TempK <- Sample.IDs.In.Bicluster[k]
            Temp.Intersect <- intersect(Sample.ID.TempK,Sample.IDs.For.GeneJ)
            if (length(Temp.Intersect) != 0){
              Genes.With.SampleK <- c(Genes.With.SampleK,Gene.IDs.In.Bicluster[j])
            }
          }
          List.Genes.With.SampleK[[k]] <- Genes.With.SampleK
        }
        
        List.With.Genes.Membership.Per.Sample.Per.Bicluster[[i]] <- List.Genes.With.SampleK
        
      }
      
      List.p.val.biclusters <- vector(mode = "list",length = nrow(Bicluster.Sample.Matrix))
      List.With.SampleIDs.Per.Bicluster <- vector(mode = "list",length = nrow(Bicluster.Sample.Matrix))
      List.With.Ser.Nos.Per.Bicluster <- vector(mode = "list",length = nrow(Bicluster.Sample.Matrix))
      for (i in 1:length(List.With.Genes.Membership.Per.Sample.Per.Bicluster))
      {
        Temp.ListI <- List.With.Genes.Membership.Per.Sample.Per.Bicluster[[i]]
        Gene.IDs.In.Bicluster <- Gene.Degrees.Bicluster.df$Gene.ID[which(Gene.Degrees.Bicluster.df$Bicluster.No == i)]
        
        p.val.Overlaps <- vector(mode = "numeric",length = length(Temp.ListI))
        
        for (j in 1:length(Temp.ListI))
        {
          t <- length(intersect(Temp.ListI[[j]],Gene.IDs.In.Bicluster))
          n <- length(Gene.IDs.Vec)
          a <- length(Gene.IDs.In.Bicluster)
          b <- a
          
          p.val.Overlaps[j] <- sum(dhyper(t:b,a,n-a,b))
        }
        List.p.val.biclusters[[i]] <- p.val.Overlaps
        List.With.SampleIDs.Per.Bicluster[[i]] <- Sample.IDs[which(Bicluster.Sample.Matrix[i,] == 1)]
        List.With.Ser.Nos.Per.Bicluster[[i]] <- rep(i,length(List.With.SampleIDs.Per.Bicluster[[i]]))
      }
      
    } else { 
      List.With.Genes.Membership.Per.Sample.Per.Bicluster <- vector(mode = "list",length = nrow(Bicluster.Sample.Matrix))
      for (i in 1:nrow(Bicluster.Sample.Matrix))
      {
        Ser.Nos.Samples.In.Bicluster <- which(Bicluster.Sample.Matrix[i,] == 1)
      
        Sample.IDs.In.Bicluster <- Sample.IDs[Ser.Nos.Samples.In.Bicluster]
      
        Start.Index <- ncol(Exp.Matrix) - length(Ser.Nos.Samples.In.Bicluster) +1
      
        Gene.IDs.In.Bicluster <- Gene.Degrees.Bicluster.df$Gene.ID[which(Gene.Degrees.Bicluster.df$Bicluster.No == i)]
      
        List.Genes.With.SampleK <- vector(mode = "list",length = length(Sample.IDs.In.Bicluster))
        for (k in 1:length(Sample.IDs.In.Bicluster))
        {
          Genes.With.SampleK <- vector(mode = "character")
          for (j in 1:length(Gene.IDs.In.Bicluster))
          {
            if (length(which(Gene.IDs.Vec == Gene.IDs.In.Bicluster[j])) > 1){
              TempJ <- which(Gene.IDs.Vec == Gene.IDs.In.Bicluster[j])[1]
            } else {
              TempJ <- which(Gene.IDs.Vec == Gene.IDs.In.Bicluster[j])
            }
          
            if (Lengths.Non.Zero.Samples[TempJ] < length(Sample.IDs.In.Bicluster)){
              Sample.IDs.For.GeneJ <- Sample.IDs[which(Exp.Matrix[TempJ,] != 0)]
            } else {
              Sample.IDs.For.GeneJ <- Sample.IDs[Sample.Order.Increasing.Matrix[TempJ,Start.Index:ncol(Exp.Matrix)]]
            }
          
            Sample.ID.TempK <- Sample.IDs.In.Bicluster[k]
            Temp.Intersect <- intersect(Sample.ID.TempK,Sample.IDs.For.GeneJ)
            if (length(Temp.Intersect) != 0){
              Genes.With.SampleK <- c(Genes.With.SampleK,Gene.IDs.In.Bicluster[j])
            }
          }
          List.Genes.With.SampleK[[k]] <- Genes.With.SampleK
        }
      
        List.With.Genes.Membership.Per.Sample.Per.Bicluster[[i]] <- List.Genes.With.SampleK
      
      }
      
      List.p.val.biclusters <- vector(mode = "list",length = nrow(Bicluster.Sample.Matrix))
      List.With.SampleIDs.Per.Bicluster <- vector(mode = "list",length = nrow(Bicluster.Sample.Matrix))
      List.With.Ser.Nos.Per.Bicluster <- vector(mode = "list",length = nrow(Bicluster.Sample.Matrix))
      for (i in 1:length(List.With.Genes.Membership.Per.Sample.Per.Bicluster))
      {
        Temp.ListI <- List.With.Genes.Membership.Per.Sample.Per.Bicluster[[i]]
        Gene.IDs.In.Bicluster <- Gene.Degrees.Bicluster.df$Gene.ID[which(Gene.Degrees.Bicluster.df$Bicluster.No == i)]
        
        p.val.Overlaps <- vector(mode = "numeric",length = length(Temp.ListI))
        
        for (j in 1:length(Temp.ListI))
        {
          t <- length(intersect(Temp.ListI[[j]],Gene.IDs.In.Bicluster))
          n <- length(Gene.IDs.Vec)
          a <- length(Gene.IDs.In.Bicluster)
          b <- a
          
          p.val.Overlaps[j] <- sum(dhyper(t:b,a,n-a,b))
        }
        List.p.val.biclusters[[i]] <- p.val.Overlaps
        List.With.SampleIDs.Per.Bicluster[[i]] <- Sample.IDs[which(Bicluster.Sample.Matrix[i,] == 1)]
        List.With.Ser.Nos.Per.Bicluster[[i]] <- rep(i,length(List.With.SampleIDs.Per.Bicluster[[i]]))
      }
      
    }
    
    Vec.p.val.biclusters <- unlist(List.p.val.biclusters)
    Vec.FDR.biclusters <- p.adjust(Vec.p.val.biclusters,method = "BH")
    
    Vec.Sample.IDs.biclusters <- unlist(List.With.SampleIDs.Per.Bicluster) 
    
    Vec.Ser.Nos.Biclusters <- unlist(List.With.Ser.Nos.Per.Bicluster)
    
    Sample.Enrichment.df <- data.frame(Vec.Ser.Nos.Biclusters,Vec.Sample.IDs.biclusters,Vec.p.val.biclusters,Vec.FDR.biclusters)
    colnames(Sample.Enrichment.df) <- c("Bicluster.No","Sample.ID","p.val","FDR(BH)")
    
    Matrix.With.P.Values <- matrix(1,nrow = max(Sample.Enrichment.df$Bicluster.No), ncol = length(Sample.IDs))
    Matrix.With.FDR.Values <- matrix(1,nrow = max(Sample.Enrichment.df$Bicluster.No), ncol = length(Sample.IDs))
    for (i in 1:length(Sample.Enrichment.df$Bicluster.No))
    {
      Matrix.With.P.Values[Sample.Enrichment.df$Bicluster.No[i],which(Sample.IDs == Sample.Enrichment.df$Sample.ID[i])] <- Sample.Enrichment.df$p.val[i]
      Matrix.With.FDR.Values[Sample.Enrichment.df$Bicluster.No[i],which(Sample.IDs == Sample.Enrichment.df$Sample.ID[i])] <- Sample.Enrichment.df$`FDR(BH)`[i]
    }
    colnames(Matrix.With.P.Values) <- Sample.IDs
    colnames(Matrix.With.FDR.Values) <- Sample.IDs
    
    pval.df <- data.frame(Bicluster.Sample.df$Bicluster.No,Matrix.With.P.Values)
    FDR.df <- data.frame(Bicluster.Sample.df$Bicluster.No,Matrix.With.FDR.Values)
    colnames(pval.df) <- c("Bicluster.No",colnames(Matrix.With.P.Values))
    colnames(FDR.df) <- c("Bicluster.No",colnames(Matrix.With.FDR.Values))
    
    FileNamePVal <- paste0(gsub(".csv","",BiclusterSamplesMatrixFileName),"_PValMatrix.csv")
    FileNameFDR <- paste0(gsub(".csv","",BiclusterSamplesMatrixFileName),"_FDRMatrix.csv")
    
    write.table(pval.df,file = FileNamePVal,row.names = F,col.names = T,sep = ",")
    write.table(FDR.df,file = FileNameFDR,row.names = F,col.names = T,sep = ",")
    
  } else if (highORlow == "low" | highORlow == "LOW" | highORlow == "l" | highORlow == "L" | highORlow == "lo" | highORlow == "LO"){
    #For LOW

    library(data.table)

    #Expression Matrix

    Exp.df <- fread(InputMatrixFileName)

    Exp.Matrix <- as.matrix(Exp.df[,-1])

    Gene.IDs.Vec <- Exp.df$Gene.ID

    #Matrix that contains the order of the samples in increasing order of expression for all genes

    Sample.Order.Increasing.Matrix <- t(apply(Exp.Matrix,1,order))

    #Read in bicluster sample file

    Bicluster.Sample.df <- fread(BiclusterSamplesMatrixFileName)

    #Matrix with sample membership in biclusters

    Bicluster.Sample.Matrix <- as.matrix(Bicluster.Sample.df[,-1])

    Sample.IDs <- colnames(Bicluster.Sample.Matrix)

    #Read in gene degrees bicluster file

    Gene.Degrees.Bicluster.df <- fread(GenesBiclustersInfoFileName)

    #Number of samples in dataset

    No.of.Samples <- ncol(Bicluster.Sample.Matrix)

    #Serial numbers of samples

    Ser.Nos.Of.Samples <- seq(1,No.of.Samples,1)

    #Ser.Nos of Non-Zero Samples per gene
    List.of.Ser.Nos.of.Non.Zero.Samples <- vector(mode = "list",length = length(Gene.IDs.Vec))
    for (i in 1:length(Gene.IDs.Vec))
    {
      List.of.Ser.Nos.of.Non.Zero.Samples[[i]] <- as.numeric(which(Exp.Matrix[i,] !=0))
    }

    Lengths.Non.Zero.Samples <- unlist(lapply(List.of.Ser.Nos.of.Non.Zero.Samples,length))

    Lengths.Zero.Samples <- No.of.Samples - Lengths.Non.Zero.Samples

    #Code to find enrichments in bottom X samples for each gene in biclusters

    if (length(which(Lengths.Non.Zero.Samples < No.of.Samples)) != 0){
      List.Of.Enrichments.Per.Gene.Per.Bicluster <- vector(mode = "list",length = nrow(Bicluster.Sample.Matrix))
      for (i in 1:nrow(Bicluster.Sample.Matrix))
      {
        Ser.Nos.Samples.In.Bicluster <- which(Bicluster.Sample.Matrix[i,] == 1)
    
        Sample.IDs.In.Bicluster <- Sample.IDs[Ser.Nos.Samples.In.Bicluster]
    
        Start.Index <- 1
    
        Gene.IDs.In.Bicluster <- Gene.Degrees.Bicluster.df$Gene.ID[which(Gene.Degrees.Bicluster.df$Bicluster.No == i)]
        List.Containing.Sample.Order <- vector(mode = "list",length = length(Gene.IDs.In.Bicluster))
    
        p.val.Overlaps <- vector(mode = "numeric",length = length(Gene.IDs.In.Bicluster))
        for (j in 1:length(Gene.IDs.In.Bicluster))
        {
          TempJ <- which(Gene.IDs.Vec == Gene.IDs.In.Bicluster[j])
          Sample.IDs.For.GeneJ <- Sample.IDs[Sample.Order.Increasing.Matrix[TempJ,Start.Index:length(Ser.Nos.Samples.In.Bicluster)]]
      
          t <- length(intersect(Sample.IDs.For.GeneJ,Sample.IDs.In.Bicluster))
          n <- length(Sample.IDs)
          a <- length(Sample.IDs.In.Bicluster)
      
          if (Lengths.Zero.Samples[TempJ] > length(Sample.IDs.In.Bicluster)){
            b <- a
          } else {
            b <- length(Sample.IDs.For.GeneJ)
          }
      
          p.val.Overlaps[j] <- sum(dhyper(t:b,a,n-a,b))
        }
        List.Of.Enrichments.Per.Gene.Per.Bicluster[[i]] <- p.val.Overlaps 
      }
    } else {
      List.Of.Enrichments.Per.Gene.Per.Bicluster <- vector(mode = "list",length = nrow(Bicluster.Sample.Matrix))
      for (i in 1:nrow(Bicluster.Sample.Matrix))
      {
        Ser.Nos.Samples.In.Bicluster <- which(Bicluster.Sample.Matrix[i,] == 1)
    
        Sample.IDs.In.Bicluster <- Sample.IDs[Ser.Nos.Samples.In.Bicluster]
    
        Start.Index <- 1
    
        Gene.IDs.In.Bicluster <- Gene.Degrees.Bicluster.df$Gene.ID[which(Gene.Degrees.Bicluster.df$Bicluster.No == i)]
        List.Containing.Sample.Order <- vector(mode = "list",length = length(Gene.IDs.In.Bicluster))
    
        p.val.Overlaps <- vector(mode = "numeric",length = length(Gene.IDs.In.Bicluster))
        for (j in 1:length(Gene.IDs.In.Bicluster))
        {
          TempJ <- which(Gene.IDs.Vec == Gene.IDs.In.Bicluster[j])
          Sample.IDs.For.GeneJ <- Sample.IDs[Sample.Order.Increasing.Matrix[TempJ,Start.Index:length(Ser.Nos.Samples.In.Bicluster)]]
      
          t <- length(intersect(Sample.IDs.For.GeneJ,Sample.IDs.In.Bicluster))
          n <- length(Sample.IDs)
          a <- length(Sample.IDs.In.Bicluster)
          b <- length(Sample.IDs.For.GeneJ)
      
          p.val.Overlaps[j] <- sum(dhyper(t:b,a,n-a,b))
        }
        List.Of.Enrichments.Per.Gene.Per.Bicluster[[i]] <- p.val.Overlaps 
      }
    }

    Vec.With.Enrichments.From.List <- unlist(List.Of.Enrichments.Per.Gene.Per.Bicluster)

    FDR.Vec.With.Enrichments.From.List <- p.adjust(Vec.With.Enrichments.From.List,method = "BH")

    New.Gene.Degrees.Bicluster.df <- data.frame(Gene.Degrees.Bicluster.df,Vec.With.Enrichments.From.List,FDR.Vec.With.Enrichments.From.List)
    colnames(New.Gene.Degrees.Bicluster.df) <- c(colnames(Gene.Degrees.Bicluster.df),"p.val.Enrichment.of.Bottom.Samples.In.Bicluster","FDR.Enrichment.of.Bottom.Samples.In.Bicluster")
    
    FileName <- paste0(gsub(".csv","",GenesBiclustersInfoFileName),"_WithAssociationPValue.csv")

    write.table(New.Gene.Degrees.Bicluster.df,file = FileName,row.names = F,col.names = T,sep = ",")
    
    #Code to find enrichments in bottom X samples for each gene in biclusters
    
    if (length(which(Lengths.Non.Zero.Samples < No.of.Samples)) != 0){
      List.With.Genes.Membership.Per.Sample.Per.Bicluster <- vector(mode = "list",length = nrow(Bicluster.Sample.Matrix))
      for (i in 1:nrow(Bicluster.Sample.Matrix))
      {
        Ser.Nos.Samples.In.Bicluster <- which(Bicluster.Sample.Matrix[i,] == 1)
        
        Sample.IDs.In.Bicluster <- Sample.IDs[Ser.Nos.Samples.In.Bicluster]
        
        Start.Index <- 1
        
        Gene.IDs.In.Bicluster <- Gene.Degrees.Bicluster.df$Gene.ID[which(Gene.Degrees.Bicluster.df$Bicluster.No == i)]
        
        
        List.Genes.With.SampleK <- vector(mode = "list",length = length(Sample.IDs.In.Bicluster))
        for (k in 1:length(Sample.IDs.In.Bicluster))
        {
          Genes.With.SampleK <- vector(mode = "character")
          for (j in 1:length(Gene.IDs.In.Bicluster))
          {
            if (length(which(Gene.IDs.Vec == Gene.IDs.In.Bicluster[j])) > 1){
              TempJ <- which(Gene.IDs.Vec == Gene.IDs.In.Bicluster[j])[2]
            } else {
              TempJ <- which(Gene.IDs.Vec == Gene.IDs.In.Bicluster[j])
            }
            
            if (Lengths.Zero.Samples[TempJ] > length(Sample.IDs.In.Bicluster)){
              Sample.IDs.For.GeneJ <- Sample.IDs[which(Exp.Matrix[TempJ,] == 0)]
            } else {
              Sample.IDs.For.GeneJ <- Sample.IDs[Sample.Order.Increasing.Matrix[TempJ,Start.Index:length(Ser.Nos.Samples.In.Bicluster)]]
            }
            
            Sample.ID.TempK <- Sample.IDs.In.Bicluster[k]
            Temp.Intersect <- intersect(Sample.ID.TempK,Sample.IDs.For.GeneJ)
            if (length(Temp.Intersect) != 0){
              Genes.With.SampleK <- c(Genes.With.SampleK,Gene.IDs.In.Bicluster[j])
            }
          }
          List.Genes.With.SampleK[[k]] <- Genes.With.SampleK
        }
        
        List.With.Genes.Membership.Per.Sample.Per.Bicluster[[i]] <- List.Genes.With.SampleK
        
      }
    } else {
      List.With.Genes.Membership.Per.Sample.Per.Bicluster <- vector(mode = "list",length = nrow(Bicluster.Sample.Matrix))
      for (i in 1:nrow(Bicluster.Sample.Matrix))
      {
        Ser.Nos.Samples.In.Bicluster <- which(Bicluster.Sample.Matrix[i,] == 1)
        
        Sample.IDs.In.Bicluster <- Sample.IDs[Ser.Nos.Samples.In.Bicluster]
        
        Start.Index <- 1
        
        Gene.IDs.In.Bicluster <- Gene.Degrees.Bicluster.df$Gene.ID[which(Gene.Degrees.Bicluster.df$Bicluster.No == i)]
        
        List.Genes.With.SampleK <- vector(mode = "list",length = length(Sample.IDs.In.Bicluster))
        for (k in 1:length(Sample.IDs.In.Bicluster))
        {
          Genes.With.SampleK <- vector(mode = "character")
          for (j in 1:length(Gene.IDs.In.Bicluster))
          {
            if (length(which(Gene.IDs.Vec == Gene.IDs.In.Bicluster[j])) > 1){
              TempJ <- which(Gene.IDs.Vec == Gene.IDs.In.Bicluster[j])[1]
            } else {
              TempJ <- which(Gene.IDs.Vec == Gene.IDs.In.Bicluster[j])
            }
            
            if (Lengths.Non.Zero.Samples[TempJ] < length(Sample.IDs.In.Bicluster)){
              Sample.IDs.For.GeneJ <- Sample.IDs[which(Exp.Matrix[TempJ,] != 0)]
            } else {
              Sample.IDs.For.GeneJ <- Sample.IDs[Sample.Order.Increasing.Matrix[TempJ,Start.Index:length(Ser.Nos.Samples.In.Bicluster)]]
            }
            
            Sample.ID.TempK <- Sample.IDs.In.Bicluster[k]
            Temp.Intersect <- intersect(Sample.ID.TempK,Sample.IDs.For.GeneJ)
            if (length(Temp.Intersect) != 0){
              Genes.With.SampleK <- c(Genes.With.SampleK,Gene.IDs.In.Bicluster[j])
            }
          }
          List.Genes.With.SampleK[[k]] <- Genes.With.SampleK
        }
        
        List.With.Genes.Membership.Per.Sample.Per.Bicluster[[i]] <- List.Genes.With.SampleK
        
      }
    }
    
    List.p.val.biclusters <- vector(mode = "list",length = nrow(Bicluster.Sample.Matrix))
    List.With.SampleIDs.Per.Bicluster <- vector(mode = "list",length = nrow(Bicluster.Sample.Matrix))
    List.With.Ser.Nos.Per.Bicluster <- vector(mode = "list",length = nrow(Bicluster.Sample.Matrix))
    for (i in 1:length(List.With.Genes.Membership.Per.Sample.Per.Bicluster))
    {
      Temp.ListI <- List.With.Genes.Membership.Per.Sample.Per.Bicluster[[i]]
      Gene.IDs.In.Bicluster <- Gene.Degrees.Bicluster.df$Gene.ID[which(Gene.Degrees.Bicluster.df$Bicluster.No == i)]
      
      p.val.Overlaps <- vector(mode = "numeric",length = length(Temp.ListI))
      
      for (j in 1:length(Temp.ListI))
      {
        t <- length(intersect(Temp.ListI[[j]],Gene.IDs.In.Bicluster))
        n <- length(Gene.IDs.Vec)
        a <- length(Gene.IDs.In.Bicluster)
        b <- a
        
        p.val.Overlaps[j] <- sum(dhyper(t:b,a,n-a,b))
      }
      List.p.val.biclusters[[i]] <- p.val.Overlaps
      List.With.SampleIDs.Per.Bicluster[[i]] <- Sample.IDs[which(Bicluster.Sample.Matrix[i,] == 1)]
      List.With.Ser.Nos.Per.Bicluster[[i]] <- rep(i,length(List.With.SampleIDs.Per.Bicluster[[i]]))
    }
    
    Vec.p.val.biclusters <- unlist(List.p.val.biclusters)
    Vec.FDR.biclusters <- p.adjust(Vec.p.val.biclusters,method = "BH")
    
    Vec.Sample.IDs.biclusters <- unlist(List.With.SampleIDs.Per.Bicluster) 
    
    Vec.Ser.Nos.Biclusters <- unlist(List.With.Ser.Nos.Per.Bicluster)
    
    Sample.Enrichment.df <- data.frame(Vec.Ser.Nos.Biclusters,Vec.Sample.IDs.biclusters,Vec.p.val.biclusters,Vec.FDR.biclusters)
    colnames(Sample.Enrichment.df) <- c("Bicluster.No","Sample.ID","p.val","FDR(BH)")
    
    Matrix.With.P.Values <- matrix(1,nrow = max(Sample.Enrichment.df$Bicluster.No), ncol = length(Sample.IDs))
    Matrix.With.FDR.Values <- matrix(1,nrow = max(Sample.Enrichment.df$Bicluster.No), ncol = length(Sample.IDs))
    for (i in 1:length(Sample.Enrichment.df$Bicluster.No))
    {
      Matrix.With.P.Values[Sample.Enrichment.df$Bicluster.No[i],which(Sample.IDs == Sample.Enrichment.df$Sample.ID[i])] <- Sample.Enrichment.df$p.val[i]
      Matrix.With.FDR.Values[Sample.Enrichment.df$Bicluster.No[i],which(Sample.IDs == Sample.Enrichment.df$Sample.ID[i])] <- Sample.Enrichment.df$`FDR(BH)`[i]
    }
    colnames(Matrix.With.P.Values) <- Sample.IDs
    colnames(Matrix.With.FDR.Values) <- Sample.IDs
    
    pval.df <- data.frame(Bicluster.Sample.df$Bicluster.No,Matrix.With.P.Values)
    FDR.df <- data.frame(Bicluster.Sample.df$Bicluster.No,Matrix.With.FDR.Values)
    colnames(pval.df) <- c("Bicluster.No",colnames(Matrix.With.P.Values))
    colnames(FDR.df) <- c("Bicluster.No",colnames(Matrix.With.FDR.Values))
    
    FileNamePVal <- paste0(gsub(".csv","",BiclusterSamplesMatrixFileName),"_PValMatrix.csv")
    FileNameFDR <- paste0(gsub(".csv","",BiclusterSamplesMatrixFileName),"_FDRMatrix.csv")
    
    write.table(pval.df,file = FileNamePVal,row.names = F,col.names = T,sep = ",")
    write.table(FDR.df,file = FileNameFDR,row.names = F,col.names = T,sep = ",")
  }
}
