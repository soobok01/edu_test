install.packages("ggplot2") 
library(ggplot2)

# 2) PCA 수행 ---------------------------------------------------------------
# train_X_sel: 샘플×변이 matrix or data.frame, train_y: 길이 n 샘플별 인종 레이블
pca_res <- prcomp(train_X_sel,        # 데이터
                  center = TRUE,      # 변수별 평균 0으로 중앙화
                  scale. = TRUE)      # 분산이 1이 되도록 스케일링

# 3) PCA 결과를 데이터프레임으로 정리 ----------------------------------------
scores <- as.data.frame(pca_res$x)    # prcomp 결과의 scores (샘플×PC)
scores$RACE <- factor(train_y)        # 인종 레이블을 마지막 열로 추가

# 4) PC1 vs PC2 산점도 그리기 ------------------------------------------------
ggplot(scores, aes(x = PC1, y = PC2, color = RACE)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Training Set PCA (PC1 vs PC2)",
    x     = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100, 1), "%)"),
    y     = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100, 1), "%)")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  )