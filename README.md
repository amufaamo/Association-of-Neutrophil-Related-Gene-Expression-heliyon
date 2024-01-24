# Association-of-Neutrophil-Related-Gene-Expression-heliyon
## Fig.1.
```R
library(edgeR)
library(tidyverse)
# Load count dataframe
all_ <- readRDS('fig1_jamb_df.rds')
# Load metadata which is id, age, and sex 
id <- read.csv('fig1_jamb_table.csv', header=FALSE)
names(id) <- c('id', 'age', 'sex')

# make edgeR function
myedgeR <- function(dataframe, group, design){
    count <- mutate_all(dataframe, function(x) as.numeric(as.character(gsub(",", "", x))))
    count <- as.matrix(count)
    group <- factor(group)
    d <- DGEList(counts = count, group = group)    
    keep <- filterByExpr(d, group=group)
    d <- d[keep, , keep.lib.sizes=FALSE]
    d <- calcNormFactors(d) 
    d <- estimateGLMCommonDisp(d, design)
    d <- estimateGLMTrendedDisp(d, design)
    d <- estimateGLMTagwiseDisp(d, design)
    return(d)
}
                        
# make group furyo(NWV), ryokou(WV), healthy
str_split(names(all_), pattern='_', simplify=TRUE) %>% .[,1] -> group

age <- id$age
sex <- id$sex

# make design matrix                        
design <- model.matrix(~ group + age + sex)

# normalize
myedgeR(all_, group, design) -> norm

# calculate cpm 
cpm(norm) -> all_normalized

names_group <- c('COV102',
'COV112',
'COV27',
'COV53',
'COV65',
'COV70',
'COV75',
'COV81',
'COV97',
'COV1',
'COV108',
'COV33',
'COV38',
'COV41',
'COV43',
'COV45',
'COV48',
'COV57',
'COV6',
'COV68',
'COV72',
'COV73',
'COV84',
'COV86',
'COV87',
'COV90',
'COV96',
'H1',
'H2',
'H7',
'H8',
'H9',
'H11',
'H14',
'H15',
'H19',
'H20',
'H21',
'H22',
'H23',
'H24',
'H25',
'H16')

                        
colnames(all_normalized) <- names_group
                        
# make histogram
rho <- cor(all_normalized, method="spearman")
d <- as.dist(1-rho)
h <- as.dendrogram(hclust(d, method = "ward.D"))
                        
# Fig 1
plot(h)
```

## Fig.2
```R
data <- readRDS('fig2_jamb_df.rds')
id <- read.csv('fig2_jamb_singlecell_metadata.csv', header=F)
names(id) <- c('id', 'sex', 'age')

myedgeR <- function (dataframe, group, design) 
{
    count <- mutate_all(dataframe, function(x) as.numeric(as.character(gsub(",", 
        "", x))))
    count <- as.matrix(count)
    group <- factor(group)
    d <- DGEList(counts = count, group = group)
    keep <- filterByExpr(d, group = group)
    d <- d[keep, , keep.lib.sizes = FALSE]
    d <- calcNormFactors(d)
    d <- estimateGLMCommonDisp(d, design)
    d <- estimateGLMTrendedDisp(d, design)
    d <- estimateGLMTagwiseDisp(d, design)    
    return(cpm(d))
}
group <- c('healthy', 'healthy', 'mild', 'mild', 'severe', 'severe', 'healthy','healthy','healthy','healthy', 'severe')
                        
age <- id$age
sex <- id$sex

# make design matrix                        
design <- model.matrix(~ group + age + sex)

# normalize
myedgeR(data, group, design) -> normalized

rho <- cor(normalized, method = "spearman")
d <- as.dist(1 - rho)
h <- as.dendrogram(hclust(d, method = "ward.D"))

plot(h)
```

## Fig. 3A
```R
# Download the seurat object below. I cannot upload because this file is big.
https://drive.google.com/file/d/13pLBTUr1wkpl3mjofjX3_NIAtmwom1x3/view?usp=drive_link

data <- readRDS('fig3_4_jamb_sctransform.rds')
library(Seurat)
Idents(data) <- 'condition1'
markers <- FindMarkers(data, ident.1 = "critical", ident.2 = "severe", latent.vars = c("age", "sex"), test.use = "MAST")
library(clusterProfiler)
source('Rscript/rprofile.r')
genes <- markers %>% dplyr::filter(p_val_adj < 0.1 & avg_log2FC > 0) %>% rownames()
source("Rscript/convert_geneID.r")
convert_geneID(genes, "SYMBOL", "ENTREZID") -> entre
entre$ENTREZID -> entre
ego2 <- enrichGO(gene         = entre,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
barplot(ego2)
```
## Fig4
```R
data <- readRDS('fig3_4_jamb_sctransform.rds')
DimPlot(data, group.by='seurat_clusters', label=TRUE)
DimPlot(data, group.by='condition4')


create_df <- function(data){
    Idents(data) <- 'condition4'

SplitObject(data, split.by='condition4') -> data2
markerlist <- c()
map(data2, function(x){try(FindMarkers(x, ident.1 = 'critical', ident.2 = 'severe', group.by='condition1', latent.vars = c("age", "sex"), test.use = 'MAST'))}) -> markerlist

d <- function(x){
    x <- tibble::rownames_to_column(x, 'gene')
    x <- dplyr::filter(x, !grepl('IG', gene))
    return(x)
}

map(markerlist, d) -> markerlist

lists <- c(4332,2357,6282,567,25801,6279,7305,2495,2512,3579,2215,6280,6283,978,3107,3106,5265,1520,23406,9535,2207,6402,25798,2212,6386,728,8826,5879,391,976,11031,10092,7097,6813,5788,11025,29108,10487,1604,1535,1992,387,527,3689,8635,226,1265)
source('Rscript/convert_geneID.r')
convert_geneID(lists, 'ENTREZID', 'SYMBOL')$SYMBOL -> neutrophil_related_genes

df <- data.frame()

df <- as.data.frame(neutrophil_related_genes)

modifym <- function(x){
    mutate(x, foldchange = case_when(
        p_val_adj >= 0.1　~ 0,
        p_val_adj < 0.1 ~ avg_log2FC
    ))    -> x
    x %>% dplyr::filter(gene %in% neutrophil_related_genes) -> x
    x <- x %>% dplyr::select(gene, foldchange)
    return(x)
}

map(markerlist, function(x){
    mutate(x, foldchange = case_when(
        p_val_adj >= 0.1　~ 0,
        p_val_adj < 0.1 ~ avg_log2FC
    ))    -> x
    x %>% dplyr::filter(gene %in% neutrophil_related_genes) -> x
    x %>% dplyr::select(gene, foldchange)
}) -> marker_modify2

dplyr::left_join(df, marker_modify2$'Pre-neu', by=c('neutrophil_related_genes' = 'gene')) -> df
dplyr::left_join(df, marker_modify2$'Pro-neu', by=c('neutrophil_related_genes' = 'gene')) -> df
dplyr::left_join(df, marker_modify2$'Mature', by=c('neutrophil_related_genes' = 'gene')) -> df

names(df) <- c('neutrophil_related_genes', 'Pre-neu', 'Pro-neu', 'Mature')
df <- mutate_all(df, ~replace(., is.na(.), 0))
return(df)    
}

df <- create_df(data)
```
