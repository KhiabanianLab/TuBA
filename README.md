# TuBA

New R package of TuBA - Analyzes large graphs to identify biclusters much faster than the [older version](https://github.com/Amartya101/TuBA-Older_Slower). 

## Getting Started

You will need an R version 3.4.0 or more recent in order to use these functions

### Prerequisites

You need the **data.table** package. You can install it using the following command:

```
install.packages("data.table",dependencies = T)
```

### Installing

Currently, the development version of **TuBA** can be installed using *install_github* from the **devtools** package.

Install **devtools** using:

```
install.packages("devtools",dependencies = T)
```
Once **devtools** is installed, run the following to install **TuBA**:
```
devtools::install_github("Amartya101/TuBA")
```


## Instructions for using TuBA

There are 3 functions in TuBA, which need to be employed in a sequential manner on the data set of interest. Make sure before running the first function (*DataCleaning*) that the data exists in a tabular format (either a .csv or a .txt file), in which the genes (or more generally, features) are along the rows and the samples are along the columns. The first column in the file must contain the IDs or names of the genes (no duplicates allowed). We have included an example file with RP genes in the [ToyExamples](https://github.com/Amartya101/ToyExamples) repository ("*RPGenes.csv*") for reference. 

The descriptions of the functions and instructions on how to use them are provided below:

### DataCleaning
The *DataCleaning* function requires one input - *File*. *File* should be the full name of the input file including the specifier for the format of the file (.csv or .txt). It is expected that the dataset in the input file is in the rectangular format, with genes along rows (first column should contain all the gene names) and samples along the columns (the header should contain the sample IDs). It should contain normalized counts and should be processed by the user to remove duplicate genes, NAs etc. The *DataCleaning* function evaluates if there are genes in the dataset that do not have a single non-zero expression value in any sample and removes these genes. It renames the header of the first column as Gene.ID and generates the output file in .csv format with the name of the input file appended by "*_Cleaned.csv*". It also generates a file that contains the names of all the genes that have non-zero expression levels in some of the samples. The name of this file ends with "*_GeneNames.csv*".

Here is an example of a valid function call

```
DataCleaning(File = "RPGenes.csv")
```
This should generate 2 output files with the names "*RPGenes_GeneNames.csv*" and "*RPGenes_Cleaned.csv*" in the current working directory.

### GenePairs

The *GenePairs* function has 5 arguments - *File*, *PercSetSize*, *JcdInd*, *highORlow*, *SampleEnrichment*. *File* should be the name of the cleaned file produced by the *DataCleaning* function. *PercSetSize* requires an input in terms of percentage of the total number of samples that would constitute the extremals for each gene. For example, if you wish to look at the top 10% you should specify this parameter to be 10 (and *highRlow = "h"*, more below). *JcdInd* specifies the Jaccard index cutoff that the overlap between the percentile sets must satisfy for a gene-pair to be a part of the graph. The parameter *highORlow* specifies whether we want the percentile sets for each gene to correspond to samples with the highest expression levels ("h") or lowest expression levels ("l") respectively. The *SampleEnrichment* argument is a logical input that specifies whether we wish for samples that are over-represented in percentile sets to be filtered out (it does this if set to TRUE). By default it is set to FALSE (recommended). 

The *GenePairs* function generates 2 output files. The first file contains all the gene-pairs that have significant overlaps between their percentile sets (Jaccard indices greater than *JcdInd* specified by the user). The name of this file ends with "*_GenePairs.csv*". This file does not contain the full gene IDs or gene names, instead the genes are labelled by their serial number in the input file. The second file contains a binary matrix with genes along the rows and samples along the columns. The first column contains the gene IDs or gene names obtained from the input file. For each gene (row), the presence of a sample in the percentile set is denoted by a 1, while samples not in the percentile set have 0. 

Examples of valid function calls are provided below (here we directly used *RPGenes.csv*, since it was already clean):
```
# For high expression
GenePairs(File = "RPGenes.csv",PercSetSize = 5,JcdInd = 0.2,highORlow = "h")
# For low expression
GenePairs(File = "RPGenes.csv",PercSetSize = 5,JcdInd = 0.2,highORlow = "l")
```
The first one will generate the following 2 files: "*RPGenes_H0.05_JcdInd0.2_GenePairs.csv*" and "*RPGenes_H0.05_JcdInd0.2_GenesSamples_BinaryMatrix.csv*". The notations in the middle of their names indicate the following: *H0.05* indicates that the gene-pairs were obtained for high expression ("H") with the percentile set size of 5% (0.05); *JcdInd0.2* indicates that a Jaccard index threshold of 0.2 was chosen to shortlist the gene-pairs.

The second one will generate the following 2 files: "*RPGenes_L0.05_JcdInd0.2_GenePairs.csv*" and "*RPGenes_L0.05_JcdInd0.2_GenesSamples_BinaryMatrix.csv*". The annotations in the middle of their names indicate the following: *L0.05* indicates that the gene-pairs were obtained for low expression ("L") with the percentile set size of 5% (0.05); *JcdInd0.2* indicates that a Jaccard index threshold of 0.2 was chosen to shortlist the gene-pairs.

### Biclustering

The *Biclustering* function has 5 arguments of which 2 inputs are necessary - *VariablePairs* and *BinaryMatrix*. *VariablePairs* requires the name of the gene-pairs .csv file generated by the *GenePairs* function; *BinaryMatrix* requires the name of the genes-samples binary matrix .csv file also generated by the *GenePairs* function. The 3 optional arguments include - *MinGenes*, *MinSamples*, and *SampleEnrichment*. *MinGenes* specifies the minimum number of genes desired by the user to be present in a bicluster (default is 3); *MinSamples* specifies the minimum number of samples desired by the user to be present in a bicluster (default is 2, but typically the biclusters have much more samples); *SampleEnrichment* is used to specify the level of enrichment that each sample must exhibit within a given bicluster as compared to its presence in the overall graph. The values for *SampleEnrichment* can range between 0 (excluding 0) and 1. Think of *SampleEnrichment* as a *p*-value - smaller values would require stronger sample associations with the bicluster, and would result in fewer samples in the biclusters. The default *SampleEnrichment* is 1.

*Biclustering* function generates 3 files. The file whose name ends with “*_GenesInBiclusters.csv*” contains list of genes contained in each bicluster. It also provides information about the total number of samples present in the bicluster (column 3), as well how many of these samples were present in the percentile sets of each gene within the bicluster (column 4). The file that ends with “*_BiclusterSamplesMatrix.csv*” is a binary matrix with biclusters along the rows and samples along the columns. For each bicluster, only the samples that are present in the bicluster are given a value of 1. The file ending with “*_GenesBiclusterSamplesMatrix.csv*” contains the genes present in each bicluster along with information about which samples are contributed to the bicluster by each gene (1 for samples contributed by a gene, 0 otherwise). The first column of this file contains the gene IDs/names, while the second column contains the bicluster number that these genes belong to; the rest of the file is comprised of the binary matrix.

Example of a valid function call is provided below:
```
Biclustering(VariablePairs = "RPGenes_H0.05_JcdInd0.2_GenePairs.csv",BinaryMatrix = "RPGenes_H0.05_JcdInd0.2_GenesSamples_BinaryMatrix.csv")
```
This will generate the following output files: "*RPGenes_H0.05_JcdInd0.2_MinGenesX_MinSamplesY_GenesInBiclusters.csv*", "*RPGenes_H0.05_JcdInd0.2_MinGenesX_MinSamplesY_BiclusterSamplesMatrix.csv*", and "*RPGenes_H0.05_JcdInd0.2_MinGenesX_MinSamplesY_GenesBiclusterSamplesMatrix.csv*". The additional annotations in the middle of their names indicate the following: *MinGenesX* indicates the number of genes in the bicluster with the fewest genes (for the example with default choices above, *X* will be 3), and *MinSamplesY* indicates the number of samples (*Y*) in the bicluster that has the fewest samples.

## Authors

* **Amartya Singh** - [Amartya101](https://github.com/Amartya101/)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hossein Khiabanian [KhiabanianLab](https://github.com/KhiabanianLab/)
* Gyan Bhanot
* Lodovico Terzi di Bergamo
