```{r format_extant}
# Load data into variable
extant_raw<-read.csv("~/Desktop/ESPM 108B/proj/rhino_extant.csv")
# Rename column names to make it more readable
colnames(extant_raw)[3]<-'DB1.1'
colnames(extant_raw)[5]<-'DB23.1'
colnames(extant_raw)[7]<-'DB44.1'
colnames(extant_raw)[9]<-'BR4.1'
colnames(extant_raw)[11]<-'BR6.1'
colnames(extant_raw)[13]<-'BR17.1'
colnames(extant_raw)[15]<-'SW35.1'
colnames(extant_raw)[17]<-'RHI32A.1'
colnames(extant_raw)[19]<-'DB14.1'
colnames(extant_raw)[21]<-'B1RH2B.1'
colnames(extant_raw)[23]<-'B1RH37D.1'
# Make a dataframe to store genotype data, one column for each locus
rhino_extant<-data.frame(matrix(nrow = nrow(extant_raw), ncol = 11))
# Fill dataframe with concatenated alleles
rhino_extant[1]<-as.numeric(paste(extant_raw$DB1, extant_raw$DB1.1, sep=""))
rhino_extant[2]<-as.numeric(paste(extant_raw$DB23, extant_raw$DB23.1, sep=""))
rhino_extant[3]<-as.numeric(paste(extant_raw$DB44, extant_raw$DB44.1, sep=""))
rhino_extant[4]<-as.numeric(paste(extant_raw$BR4, extant_raw$BR4.1, sep=""))
rhino_extant[5]<-as.numeric(paste(extant_raw$BR6, extant_raw$BR6.1, sep=""))
rhino_extant[6]<-as.numeric(paste(extant_raw$BR17, extant_raw$BR17.1, sep=""))
rhino_extant[7]<-as.numeric(paste(extant_raw$SW35, extant_raw$SW35.1, sep=""))
rhino_extant[8]<-as.numeric(paste(extant_raw$RHI32A, extant_raw$RHI32A.1, sep=""))
rhino_extant[9]<-as.numeric(paste(extant_raw$DB14, extant_raw$DB14.1, sep=""))
rhino_extant[10]<-as.numeric(paste(extant_raw$B1RH2B, extant_raw$B1RH2B.1, sep=""))
rhino_extant[11]<-as.numeric(paste(extant_raw$B1RH37D, extant_raw$B1RH37D.1, sep=""))
#Remove rows with all NA's
rhino_extant<-rhino_extant[rowSums(is.na(rhino_extant)) != ncol(rhino_extant),]
```


```{r PCA_extant}
# Load library for PCA analysis
library(adegenet)
# Converting dataframe to 'genind' data type
x_extant<-df2genind(rhino_extant, ncode=3, NA.char = NA, ploidy = 2)
# Quantifies number of PCA components and clusters
grp_extant<-find.clusters(x_extant, max.n.clust = 20) # 200, 10
# Plots discriminant analysis eigenvalues, represents diversity between pre-defined groups
dapc1_extant <- dapc(x_extant, grp_extant$grp) # 10, 3
1# Scatter plot of diversity between three clusters
scatter(dapc1_extant, ratio.pca=0.3, bg="white", pch=20,  cell=0,
        cstar=0, solid=.4, cex=3, clab=0,
        mstree=TRUE, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, txt.leg=paste("Cluster",1:3))
# Bar plot to represent the group assignment probability of individuals to several groups
compoplot(dapc1_extant, posi="bottomright",legend = FALSE, lab="",
          ncol=1, xlab="individuals", col=funky(5))
```