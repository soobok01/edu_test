

library(randomForest)


##setwd('~~')


train_X <- read.table("TrainingSNV_set.txt", header=TRUE, sep="\t", check.names=FALSE)
train_y <- scan("TrainingRACE_Y.txt", what=character(), sep="\n")
test_X <- read.table("TestSNV_set.txt", header=TRUE, sep="\t", check.names=FALSE)
test_y <- scan("TestRACE_Y.txt", what=character(), sep="\n")

train_mat <- as.matrix(train_X)
test_mat <- as.matrix(test_X)


pops <- unique(train_y)

# AF 행렬 계산 --------------------------------------------------
AF_mat <- sapply(pops, function(pop) {
  idx <- which(train_y == pop)                     # 해당 인종 샘플 인덱스
  # genotype 값을 ALT allele count 로 바로 쓴다고 가정 (0/1/2)
  rowSums(train_mat[, idx]) / (2 * length(idx))
})

colnames(AF_mat) <- pops
rownames(AF_mat) <- rownames(train_X)

#  각 행(variant)에 대해 5개 중 하나라도 AF>=0.6 이면 채택
keep_variants <-apply(AF_mat, 1, function(x) any(x >= 0.6))

length(keep_variants)  # 선택된 variant 개수

train_sel_mat <- train_mat[keep_variants, , drop = FALSE]
test_sel_mat  <- test_mat[ keep_variants, , drop = FALSE]


# 2) 일반적인 머신러닝 함수에 맞게 sample×feature 로 전치
train_X_sel <- t(train_sel_mat)   # 이제 행 = 샘플, 열 = 변이(feature)
test_X_sel  <- t(test_sel_mat)

# 3) data.frame으로 변환하고 레이블 붙이기
train_df <- as.data.frame(train_X_sel)
train_df$RACE <- factor(train_y, levels = pops)

test_df <- as.data.frame(test_X_sel)
test_df$RACE <- factor(test_y,  levels = pops)


# 랜덤 포레스트 모델 -----------------------------------------------------
set.seed(123)  # 재현성을 위해 시드 고정
rf_model <- randomForest(
  RACE ~ .,            # formula: RACE를 나머지 변수들로 예측
  data    = train_df,  # 학습 데이터
  ntree   = 200,       # 트리 개수
  importance = TRUE    # 변수 중요도 계산
)
print(rf_model)


rf_pred <- predict(rf_model, newdata = test_df)

# 성능 평가: Confusion Matrix & Accuracy
conf_mat_rf <- table(Predicted = rf_pred, Actual = test_df$RACE)
print(conf_mat_rf)
acc_rf <- mean(rf_pred == test_df$RACE)
cat("Random Forest Accuracy:", round(acc_rf, 3), "\n\n")

##############결과 저장장
test_pred <- predict(rf_model, newdata = test_df)
test_res  <- data.frame(
  Sample    = rownames(test_df),
  Actual    = test_df$RACE,
  Predicted = test_pred,
  stringsAsFactors = FALSE
)
write.table(
  test_res,
  file      = "rf_test_predictions.txt",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

######

# 1) 변수 중요도 추출 및 저장 -----------------------------------------------
# importance() 로 얻은 행렬을 data.frame으로 변환
imp_mat <- importance(rf_model)          
imp_df  <- as.data.frame(imp_mat)        
imp_df$Feature <- rownames(imp_df)       # 컬럼명으로 있던 feature 이름을 새 컬럼으로

# 컬럼 순서: Feature, MeanDecreaseGini, MeanDecreaseAccuracy 등
imp_df <- imp_df[, c("Feature", setdiff(names(imp_df), "Feature"))]

# 파일로 쓰기 (탭 구분, row.names 제외)
write.table(
  imp_df,
  file      = "rf_feature_importance.txt",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
