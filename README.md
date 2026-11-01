# 数据清洗排纳
library(dplyr)
#1. 载入原始数据集
data(NHANES4)
#2. 定义 CKD-EPI 2021（无种族修正版）计算函数
calc_egfr_2021 <- function(scr, age, sex) {
  egfr <- numeric(length(scr))
  for(i in seq_along(scr)) {
    s <- scr[i]
    a <- age[i]
    sx <- sex[i]
    if(is.na(s) || is.na(a) || is.na(sx)) {
      egfr[i] <- NA
      next
    }
    if(sx == 2) { # 女性
      kappa <- 0.7; alpha <- -0.241; mult <- 1.012
      egfr[i] <- 142 * (min(s/kappa, 1)^alpha) * (max(s/kappa, 1)^-1.200) * (0.9938^a) * mult
    } else { # 男性
      kappa <- 0.9; alpha <- -0.302; mult <- 1.000
      egfr[i] <- 142 * (min(s/kappa, 1)^alpha) * (max(s/kappa, 1)^-1.200) * (0.9938^a) * mult
    }
  }
  return(egfr)
}
#3. 初始化记录初始样本量
n_current <- nrow(NHANES4)
cat("=== 队列排纳过程统计 ===\n")
cat("原始总样本量 (N_total):", n_current, "\n")
#新增 egfr 变量
nhanes_work <- NHANES4 %>%
  mutate(egfr = calc_egfr_2021(scr = creat, age = age, sex = gender))
#步骤 1：年龄筛选 (纳入 >= 45 岁中老年人群)
step1_data <- nhanes_work %>% filter(!is.na(age) & age >= 45)
N1 <- nrow(step1_data)
cat("步骤 1 (筛选年龄 >= 45 岁): 剔除", n_current - n_step1, "人 | 剩余:", n_step1, "人\n")
n_current <- N1
#步骤 2：排除孕妇 (pregnant == 1)
step2_data <- step1_data %>% filter(pregnant != 1 | is.na(pregnant))
N2 <- nrow(step2_data)
cat("步骤 2 (排除孕妇): 剔除", n_current - N2, "人 | 剩余:", N2, "人\n")
n_current <- N2
#步骤 3：剔除核心原生变量缺失者 (gender, creat, egfr, phenoage0, kdm0,status,permth_int)
step3_data <- step2_data %>% filter(
  !is.na(gender),
  !is.na(creat),
  !is.na(egfr),
  !is.na(phenoage0),
  !is.na(kdm0),
  !is.na(status),
  !is.na(permth_int)
)
N3 <- nrow(step3_data)
cat("步骤 3 (剔除核心变量缺失者): 剔除", n_current - N3, "人 | 剩余:", N3, "人\n")
n_current <- N3
#步骤 4：剔除极端离群值 (eGFR < 15)
nhanes_analytic <- step3_data %>% filter(egfr >= 15)
N4 <- nrow(nhanes_analytic)
cat("步骤 4 (剔除 eGFR < 15 极端值): 剔除", n_current -N4, "人 | 最终分析样本:", N4, "人\n")
cat("========================\n")

# 提取聚类变量并剔除缺失值（使用 complete.cases 锁定精确的行）
install.packages("mclust")
library(mclust)
cluster_vars <- c("phenoage0", "kdm0", "egfr", "bmi", "sbp")
complete_idx <- complete.cases(nhanes_analytic[, cluster_vars])
#提取用于聚类的干净数据子集
cluster_data_clean <- nhanes_analytic[complete_idx, cluster_vars]
#Z-score 标准化
cluster_matrix <- scale(cluster_data_clean)
#运行 GMM 聚类
gmm_fit <- Mclust(cluster_matrix)
#查看 BIC 曲线
plot(gmm_fit, what = "BIC")
#安全地将聚类标签赋给对应的干净子集（或整张表对应的行）
nhanes_analytic_clustered <- nhanes_analytic
#先初始化为 NA，再将对应完整行的位置填入聚类标签
nhanes_analytic_clustered$cluster_gmm <- NA
nhanes_analytic_clustered$cluster_gmm[complete_idx] <- as.factor(gmm_fit$classification)
#查看各亚型人数分布
table(nhanes_analytic_clustered$cluster_gmm, useNA = "ifany")
