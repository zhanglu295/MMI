# Six scripts to calculate MMI (Multimodal Mutual Information) and MDI (Multimodal Direct Information)

## 1. dynamicprogramming.m is used for clustering the expression profiles for each gene into clusters by dynamic programming.
Parameters: 
a. Data: gene expression data. each row represents one sample, each column for each gene.
b. maxallow: the maximum cluster number allowed. The program calculates the BIC values from 1 to maxallow to choose the exact cluster number with the smallest BIC.
Values:

## 2. MMI.m is used to caluculate multimodal mutual information for gene pairs.
Parameters:
a. Both X and Y are column vectors, containing the gene expression profiles for one gene, respectively. 
b. Both clusterlabelX and clusterlabelY are column vectors, they store the clustering results for gene X and Y.  

## 3. MDIin.m is used to caluculate the covariance and inverse covariance matrix of inner MDI.
Parameters:
a. expressdata: gene expression data. each row represents one sample, each column for each gene.
b. clusterassign: the clustering result for all the genes. Each row is a sample, each column is one gene.
c. clusternum: the number of cluster for all the genes.
d. shrink: shrinkage parameter for inversion.

## 4. MDIcombine is used to calculate multimodal direct information.
Parameters:
a. data:  gene expression data. each row represents one sample, each column for each gene.
b. clusterassign: the clustering result for all the genes. Each row is a sample,
 each column is one gene.
c. clusternum: the number of cluster for all the genes.
d. MDIin_invMatrix:  inverse covariance matrix of inner MDI. The output of MDIin.m
e.  shrink: shrinkage parameter for matrix inversion.

## 5. main.m is used to test the scripts of MDI and MMI.

## 6. myfuncov.m is used to approximate the real covariace matrix by constrain optimization.



