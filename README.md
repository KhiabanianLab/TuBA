# TuBA
TuBA - *Tunable Biclustering Algorithm* for gene expression datasets is a graph-based unsupervised biclustering algorithm, customized to identify alterations in tumors. It is based on the hypothesis that gene pairs relevant to an alteration associated with dysregulated gene expression share a statistically significant number of extremal samples. TuBA has been developed in the [Khiabanian Lab](http://khiabanian-lab.org/) by Amartya Singh and Dr. Hossein Khiabanian.

## INSTRUCTIONS FOR USING TuBA’s CODE

## 1. System Requirements

It is recommended that the code for TuBA be run on a machine with at least 16GB RAM with an R version not older than 3.3.0. The run time of the code can vary between 2 - 14 hrs depending on the size of the dataset and the choice of the knobs.

## 2. Packages Required

The following packages have been utilized in TuBA’s Code: data.table, plyr and igraph.

They can be installed by using the following commands in R:

`install.packages("data.table", dependencies=TRUE)`

`install.packages("plyr”, dependencies = TRUE)`

`install.packages("igraph”, dependencies = TRUE)`

## 3. How to Use The Code

The code for TuBA consists of 3 functions:

	1. DatasetPreparationFunction
	2. SignificantGenePairsFunction
	3. BiclusteringFunction

The three have to be executed serially beginning with the DatasetPreparationFunction. The serial execution is necessary since the inputs required to execute the subsequent functions are generated by the previous function.

	1. DatasetPreparationFunction

The `DatasetPreparationFunction` requires one input - *InputFileName*. The name of the input file (.txt or .csv) should be specified in full including the specifier for the format of the file. It is expected that the dataset in the input file is in the rectangular format with genes along rows (first column should contain all the gene names) and samples along the columns (the header should contain the sample IDs). It should contain normalized counts and should be processed by the user to remove duplicate genes, NAs etc. The `DatasetPreparationFunction` evaluates if there are genes in the dataset that do not have a single non-zero expression value in any sample and eliminates these genes. It renames the header of the first column as Gene.ID and generates the output file in .csv format with the name of the input file appended by “_Cleaned.csv”. For example, if the name of the input file is “GeneExpressionDataset1.txt”, the name of the output file will be “GeneExpressionDataset1_Cleaned.csv”. It also generates a file that contains the names of all the genes that have non-zero expression values in some of the samples. The name of this file ends with “_GeneNames.csv”. 

*Eg:* `DatasetPreparationFunction(InputFileName = “GeneExpressionDataset1.txt”)`

	2. SignificantGenePairsFunction
	
The `SignificantGenePairsFunction` requires three inputs - *InputFileName*, *PercentileCutOff*, *highORlow*. *InputFileName* should be the name of the cleaned file produced by the `DatasetPreparationFunction`. *PercentileCutOff* requires an input in terms of percentage of the total sample size that should be considered for comparing the extremals for each gene. For example, if you wish the percentile set size for comparison to be 10% you should specify this parameter to be 10. The input parameter *highORlow* requires you to specify whether you want the percentile sets to correspond to high expression or low expression. Use this ONLY if the data has been obtained from an RNASeq platform. In order to specify high expression you can assign any one of these character values to *highORlow* - “h”, “H”, “high”, “High”. Similarly, for low expression you can assign any one of these character values to *highORlow* - “l”, “L”, “low”, “Low”. 

The `SignificantGenePairsFunction` generates two files. The first file contains a binary matrix with 1 corresponding to the samples in the percentile set for each gene. The name of the output file contains the name of the input file appended by “_H(*PercentileCutOff*/100)_Genes_Samples_BinaryMatrix.csv” or “_L(*PercentileCutOff*/100)_Genes_Samples_BinaryMatrix.csv”, depending on whether you specify “h” or “l”. For example, for the input file “GeneExpressionDataset1_Cleaned.csv” and a percentile cutoff of 5% for high expression (*highORlow* = “h”), the name of the output file will be “GeneExpressionDataset1_H0.05_Genes_Samples_BinaryMatrix.csv”. The second file generated by the `SignificantGenePairsFunction` contains all the gene pairs that have significant overlaps between their percentile sets (p < 0.05). The name of this file ends with “_SignificantGenePairs.csv”. The genes are labelled by their serial number in the list of genes produced by the `DatasetPreparationFunction`.

In addition to these two files, this function also generates 6 plots to help guide the choice of the significance level of overlap. The plot with its ending with “_NoOfEdges.pdf” shows the total number of edges in the graph for different levels of the significance of overlap (in terms of -log10(p)). The plot with its name ending with “_NoOfSamples.pdf” represents the total number of samples present in the graph for different levels of significance of overlaps. The plot with its name ending with “_GenesAdded.pdf” shows the the number of new genes added to the graph with each incremental drop in the level of significance of overlap. Similarly, the plot with its name ending with “_SamplesAdded.pdf” shows the number of new samples added to the graph with every incremental drop in the level of significance of overlap. The plot with its name ending with “_GenesPerEdge.pdf” shows the ratio of the number of new genes added to the graph to the number of new edges added to the graph with every incremental drop in the level of significance of overlap. Similarly, the plot with its name ending with “_SamplesPerEdge.pdf” shows the ratio of the number of new samples added to the graph to the number of new edges added to the graph with every incremental drop in the level of significance of overlap.

*Eg:* `SignificantGenePairsFunction(InputFileName = “GeneExpressionDataset1_Cleaned.csv”,PercentileCutOff = 5, highORlow = “h”)`

	3. BiclusteringFunction

The `BiclusteringFunction` requires three inputs - *FileNameGenePairs*,*FileNameBinaryMatrix*,*OverlapCutOff*. *FileNameGenePairs* requires the name of the “_SignificantGenePairs.csv” file that is produced by the SignificantGenePairsFunction. *FileNameBinaryMatrix* requires the name of the binary matrix (“_Genes_Samples_BinaryMatrix.csv”) produced by the `SignificantGenePairsFunction`. The third input - *OverlapCutOff* - requires the specification of the minimum level of significance necessary for a gene pair to be represented in our graph. This is done by specifying the modulus of the exponent of 10 for any given significance level of choice. For example, if the significance level for overlap is desired to be 10^(-10), we would specify *OverlapCutOff* = 10. If instead the desired level was 10^{-12.5}, we would specify *OverlapCutOff* = 12.5

`BiclusteringFunction` generates six files. The file whose name ends with “_Edgelist_Seeds.csv” contains the gene pair associations for every seed. The associations are represented in the form of a table with genes in column 1 associated with the corresponding genes (the ones to their immediate right) in column 2. The file whose name ends with “_SeedsSamplesMatrix.csv” is a binary matrix with seeds along the rows and samples along the columns. For each seed, only the samples that are present in the seed are given a value of 1. The file ending with “_GeneNames_Seeds.csv” contains the names of genes in each seed. This can be used for functional annotation and clustering analysis. The file that ends with “_BiclusterSamplesMatrix.csv” is a binary matrix with biclusters along the rows and samples along the columns. For each bicluster, only the samples that are present in the bicluster are given a value of 1. The file ending with “_GenesDegreesInBiclusters.csv” contains the genes present in each bicluster along with their respective degrees in the biclusters. Finally, the file ending with “_Edgelist_Biclusters.csv” contains the gene pair associations that constitute every bicluster. The associations are represented in the form of a table with genes in column 1 associated with the corresponding genes (the ones to their immediate right) in column 2. Note that the names of all the six files also contain the choice of the overlap cutoff in following format - “OverlapCutOff1E-XX”, where XX is the choice of the *OverlapCutOff* specified in the input. 

Eg: `BiclusteringFunction(FileNameGenePairs = “GeneExpressionDataset1_H0.05_SignificantGenePairs.csv”,FileNameBinaryMatrix = “GeneExpressionDataset1_H0.05_Genes_Samples_BinaryMatrix.csv”, OverlapCutOff = 15)`

