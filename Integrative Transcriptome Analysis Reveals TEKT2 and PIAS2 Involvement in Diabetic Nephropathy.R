##Integrative Transcriptome Analysis Reveals TEKT2 and PIAS2 Involvement in Diabetic Nephropathy
##Array data analysis
#Load CEL data
setwd("~/mydata/cel_data/CTL Glom")
data <- ReadAffy()
librsampleNames(data)
N=length(data)
N = length(data)
N
#QC
image(data[,1])
image(data[,2])
image(data[,3])
image(data[,4])
image(data[,5])
library(affyPLM)
library(affy)
library(RColorBrewer)
color <- brewer.pal(5,"Set3")
data.ege <-AffyRNAdeg(data)
Pset<-fitPLM(data)
Mbox(Pset,col=color,main="RLE",las=1)
boxplot(Pset,col=color,main="RLE",las=1)
plotAffyRNAdeg(data.ege,col=colors())
legend("topleft",sampleNames(data),col=color,lwd = 0.5,inset = 0.05,cex = 0.3)
eset.rma<-rma(data)
ctl_glom_exprs<-exprs(eset.rma)
probeid <-rownames(ctl_glom_exprs)
ctl_glom_exprs<-cbind(probeid,ctl_glom_exprs)
write.table(ctl_glom_exprs,file = "02__DN_GLOM_exprs.txt",sep = "\t",quote = F,row.names = F)
#sesson change
normal_kidney_expr1<-read.table("02__ctl_glom_exprs.txt",sep = "\t",header = T)
diabetic_nephropathy_expr1<- read.table("02__DN_GLOM_exprs.txt",sep = "\t",header = T)
glom_exprs<- merge(normal_kidney_expr1,diabetic_nephropathy_expr1,by ="probeid")
write.table(glom_exprs,file = "03GLOM_exprs571.txt",sep = "\t",quote = F,row.names = F)
probe_exp<-read.table("03GLOM_exprs571.txt",sep = "\t",header =T,row.names = 1)
probeid_geneid<-read.table("GPL571-17391.txt",sep = "\t",header =T)
probe_name<-rownames(probe_exp)
loc<-match(probeid_geneid[,1],probe_name)
probe_exp<-probe_exp[loc,]
raw_geneid<-as.numeric(as.matrix(probeid_geneid[,3]))
# Create a probe id index with gene id
index<-which(!is.na(raw_geneid))
# Extract the probe with the gene id
geneid<-raw_geneid[index]
# Find expression values for each gene id
exp_matrix<-probe_exp[index,]
geneidfactor<-factor(geneid)
# When multiple probes correspond to a gene, calculate the average expression level
gene_exp_matrix<-apply(exp_matrix,2, function(x) tapply(x, geneidfactor, mean))
#Use geneid as rownames
rownames(gene_exp_matrix)<-levels(geneidfactor)
geneid<-rownames(gene_exp_matrix)
gene_exp_matrix2<-cbind(geneid,gene_exp_matrix)
write.table(gene_exp_matrix2,file="10GLOM571.geneid.expr.txt",sep = "\t",quote = F,row.names = F)
# Convert geneid to genesymbol
loc<-match(rownames(gene_exp_matrix),probeid_geneid[,3])
rownames(gene_exp_matrix)=probeid_geneid[loc,2]
genesymbol<- rownames(gene_exp_matrix)
gene_exp_matrix3<- cbind(genesymbol,gene_exp_matrix)
write.table(gene_exp_matrix3,file = "11GLOM_exprs571_genesymbol.txt",sep = "\t",quote = F,row.names = F)
# imputation
library(impute)
gene_exp_matrix<- read.table("11GLOM_exprs571_genesymbol.txt",sep = "\t",header =T,row.names = 1)
gene_exp_matrix<-as.matrix(gene_exp_matrix)
imputed_gene_exp<- impute.knn(gene_exp_matrix,k=10,rowmax = 0.5,colmax = 0.8,maxp = 3000,rng.seed = 362436069)
GeneExp<-imputed_gene_exp$data
#Save the result as a table
genesymbol<- rownames(GeneExp)
GeneExp<- cbind(genesymbol,GeneExp)
write.table(GeneExp,file= "12GLOM_exprs571_genesymbol.txt",sep = "\t",quote = F,row.names = F)

##PCA analysis
A<-read.csv(file = "glom_571_PCA.csv",header = T,row.names = 1)
iris.pca <- PCA(A, graph = FALSE)
group<-read.csv(file = "glom_group_GPL571.csv",header = T)
fviz_pca_ind(iris.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = group$group, # color by groups
             palette = "jco", 
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups")
## Differential gene analysis with limma
library(limma)
eset<-read.table("12GLOM_exprs571_genesymbol.txt",sep = "\t",header =T,row.names = 1)
# Set up the number of samples in each group
condition= factor(c(rep("C",13),rep("N",9)))
#Creat matrix
design<- model.matrix(~-1+condition)
colnames(design)<- c("CTL","DN")
contranst.matrix<-makeContrasts(contrasts = "DN-CTL",levels = design)
fit<-lmFit(eset,design)
fit1<- contrasts.fit(fit,contranst.matrix)
fit2<-eBayes(fit1)
dif<- topTable(fit2,coef = "DN-CTL",n= nrow(fit2),adjust="BH")
genesymbol<- rownames(dif)
dif<- cbind(genesymbol,dif)
write.table(dif,file="14GLOM_dif_FCgene_571.txt",sep = "\t",quote = F,row.names = F)
#Find differentially expressed genes with FC greater than 1.5
dif2<- topTable(fit2,coef = "DN-CTL",n= nrow(fit2),lfc = log2(1.5),adjust="BH")
dif2<-dif2[dif2[,"P.Value"]<0.05,]
genesymbol<- rownames(dif2)
dif2<-cbind(genesymbol,dif2)

write.table(dif2,file="15__GLOM571_FC1.5FDR0.05.txt",sep = "\t",quote = F,row.names = F)

##KEGG/GO enrichment
library(org.Hs.eg.db)
library(clusterProfiler)
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
A<-read.table(file = "15__GLOM1009_FC1.5FDR0.05.txt ", sep=”\t”)
gene<-A$genesymbol
gene.df <- bitr(gene, fromType = " SYMBOL ",
                toType = c("ENSEMBL", " ENTREZID "),
                OrgDb = org.Hs.eg.db)

genelist<-gene.df$ENTREZID
#KEGG analyse
kk <- enrichKEGG(gene = genelist,
                 organism ="hsa",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 minGSSize = 1,
                 #readable = TRUE ,
                 use_internal_data =FALSE)
dotplot(kk, showCategory=20,title="EnrichmentKEGG")
#GO analyse
ego_all<- enrichGO(gene = gene.df$ENTREZID,
                   OrgDb= org.Hs.eg.db,
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
dotplot(ego_all, showCategory=20,title="EnrichmentGO")

##GSEA analysis
Follow the instruction from https://www.gsea-msigdb.org/gsea/index.jsp

##Deconvolution
Follow the instruction of https://github.com/cran/BisqueRNA







