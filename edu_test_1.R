getwd()
## set working directory 
setwd("C:/Users/KOBIC/Desktop/coad_gene_expression")


#library install
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biobase")
if (!require("NMF", quietly = TRUE))
  install.packages("NMF")

#library load
library(NMF)
library(Biobase)

# 1) data load & log2 transform
raw_data <- read.table('coad_gene_expression.txt', sep = "\t", header = TRUE, row.names = 1)
expr_log2 <- log2(raw_data + 2)

# 2) 상위 1000개 추출

sd_vals <- apply(expr_log2, 1, sd, na.rm = TRUE)
top1000_idx <- order(sd_vals, decreasing = TRUE)[1:1000]
expr_log2_top1000 <- expr_log2[top1000_idx, ]


# 3) nmf test

nmf_res_2 <- nmf(expr_log2_top1000, rank = 2, nrun = 100, seed = 123)
nmf_res_3 <- nmf(expr_log2_top1000, rank = 3, nrun = 100, seed = 123)
nmf_res_4 <- nmf(expr_log2_top1000, rank = 4, nrun = 100, seed = 123)


consensusmap(nmf_res_2)
consensusmap(nmf_res_3)
consensusmap(nmf_res_4)

layout(cbind(1,2)) 
basismap(nmf_res_3, subsetRow=TRUE) 
coefmap(nmf_res_3)

# 4)feature extraction

W <- basis(nmf_res_3)    # m × r
H <- coef(nmf_res_3)     # r × n

s<-extractFeatures(nmf_res_3)
metagene_1_idx=s[[1]]
metagene_2_idx=s[[2]]
metagene_3_idx=s[[3]]

c1_metagenes=rownames(expr_log2_top1000)[metagene_1_idx]
c2_metagenes=rownames(expr_log2_top1000)[metagene_2_idx]
c3_metagenes=rownames(expr_log2_top1000)[metagene_3_idx]

c1_metagenes <- strsplit(as.character(c1_metagenes), "\\|")
c1_metagenes_new <- do.call(rbind, c1_metagenes)
c2_metagenes <- strsplit(as.character(c2_metagenes), "\\|")
c2_metagenes_new <- do.call(rbind, c2_metagenes)
c3_metagenes <- strsplit(as.character(c3_metagenes), "\\|")
c3_metagenes_new <- do.call(rbind, c3_metagenes)


write.table(c1_metagenes_new[,1],'c1_metagene.txt',
            sep        = "\t",
            quote      = FALSE,
            row.names  = FALSE,
            col.names  = FALSE)

write.table(c2_metagenes_new[,1],'c2_metagene.txt',
            sep        = "\t",
            quote      = FALSE,
            row.names  = FALSE,
            col.names  = FALSE)

write.table(c3_metagenes_new[,1],'c3_metagene.txt',
            sep        = "\t",
            quote      = FALSE,
            row.names  = FALSE,
            col.names  = FALSE)




# 6) 샘플 클러스터 레이블 추출 (Hard clustering)
clusters <- apply(H, 2, which.max)
print(table(clusters))

write.table(clusters,'clusters.txt',
           sep        = "\t",
           quote      = FALSE,
           row.names  = TRUE,
           col.names  = FALSE)
# 7) PCA를 통한 클러스터 시각화
pca_out <- prcomp(t(expr_log2_top1000), scale. = TRUE)
plot(pca_out$x[,1:2], col = clusters, pch = 16,
     xlab = "PC1", ylab = "PC2",
     main = "PCA of Samples colored by NMF Clusters")



#############################




estim.r<-nmf(expr_log2_top1000,2:6,nrun=10,seed=123456) 

plot(estim.r)

consensusmap(estim.r)


consensusmap(estim.r[[3]])

k1=estim.r$fit$`2`


names(estim.r$fit)
