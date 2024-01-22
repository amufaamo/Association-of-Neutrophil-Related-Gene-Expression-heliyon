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

まずはGoogle documentに書いていく
https://docs.google.com/document/d/1TGLvam3WYXM4m24sBs_iaXSF9z9YlaodvyV0FA1RKZ8/edit
