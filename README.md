library(dplyr)
library(mclust)
library(tableone)
# 数据排纳
#1. 载入原始数据集
data(NHANES4)
#2. 定义 CKD-EPI 2021 计算函数
calc_egfr_2021 <- function(scr, age, sex) {
  egfr <- numeric(length(scr))
  for(i in seq_along(scr)) {
    s <- scr[i]; a <- age[i]; sx <- sex[i]
    if(is.na(s) || is.na(a) || is.na(sx)) { egfr[i] <- NA; next }
    if(sx == 2) {
      egfr[i] <- 142 * (min(s/0.7, 1)^-0.241) * (max(s/0.7, 1)^-1.200) * (0.9938^a) * 1.012
    } else {
      egfr[i] <- 142 * (min(s/0.9, 1)^-0.302) * (max(s/0.9, 1)^-1.200) * (0.9938^a) * 1.000
    }
  }
  return(egfr)
}
#3.【关键步骤】：在排纳前，利用临床标准重新衍生高血压和糖尿病状态
nhanes_work <- NHANES4 %>%
  mutate(
    egfr = calc_egfr_2021(scr = creat, age = age, sex = gender),
        # 【重新定义高血压】：平均 SBP >= 140 或 DBP >= 90，或者原本数据中自报/确诊（若有记录）
    # 这里主要基于生理测量值判定，可完美避免 93% 的高缺失率
    new_hyperten = if_else(sbp >= 140 | dbp >= 90 | hyperten == 1, 1, 0, missing = 0),
        # 【重新定义糖尿病】：HbA1c > 6.5% 或 空腹血糖 >= 7.0 mmol/L，或者原本数据标记为 1
    new_diabetes = if_else(hba1c > 6.5 | glucose_fasting >= 7.0 | diabetes == 1, 1, 0, missing = 0)
  )
#4. 规范排纳流程
nhanes_analytic <- nhanes_work %>%
  filter(!is.na(age) & age >= 45) %>%
  filter(pregnant != 1 | is.na(pregnant)) %>%
  filter(
    !is.na(gender), !is.na(creat), !is.na(egfr),
    !is.na(status), !is.na(permth_int)
  ) %>%
  filter(egfr >= 15) %>%
  filter(
    !is.na(phenoage0),
    !is.na(kdm0),
    !is.na(bmi),
    !is.na(sbp)
  )
cat("重新定义并清洗后的最终分析样本量:", nrow(nhanes_analytic), "\n")

# GMM聚类
library(dplyr)
library(mclust)
library(ggplot2) 
#前提假设：包含了 phenoage0, kdm0, egfr, bmi, sbp 以及随访变量，且无任何 NA 缺失值。
#1. 提取聚类变量并进行 Z-score 标准化（消除量纲与数值大小差异）
cluster_vars <- c("phenoage0", "kdm0", "egfr", "bmi", "sbp")
cluster_matrix <- scale(nhanes_analytic[, cluster_vars])
#2. 【全局评估阶段】计算不同群数 (G) 和协方差模型下的 BIC 值
bic_values <- mclustBIC(cluster_matrix)
#3. 【可视化阶段】使用 ggplot2 绘制高颜值 BIC 曲线图，观察拐点
bic_df <- data.frame(
  G = rep(1:nrow(bic_values), times = ncol(bic_values)),
  BIC = as.vector(bic_values),
  Model = rep(colnames(bic_values), each = nrow(bic_values))
) %>% na.omit()
ggplot(bic_df, aes(x = G, y = BIC, color = Model, group = Model)) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  geom_point(size = 1.5) +
  theme_minimal(base_size = 13) +
  labs(
    title = "GMM Model Selection via BIC (with KDM-BA)",
    subtitle = "Higher BIC value indicates better model fit",
    x = "Number of Components (Clusters, G)",
    y = "BIC Value",
    color = "Covariance Model"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = unique(bic_df$G))
  gmm_fit$BIC
#4. 【最终锁定阶段】根据 BIC 曲线拐点（确定选择 G = 4）
#【核心要点】：设置随机种子，彻底锁定 EM 算法的迭代初值，确保每次运行结果绝对一致、可复现
set.seed(123) 
optimal_G <- 4 
#运行最终的 GMM 模型
gmm_fit_final <- Mclust(cluster_matrix, G = optimal_G)
#查看最终选定模型的详细统计摘要
cat("=== 最终拟合模型摘要 ===\n")
print(summary(gmm_fit_final))
#5. 将聚类标签安全地赋回分析数据框
nhanes_analytic_clustered <- nhanes_analytic
nhanes_analytic_clustered$cluster_gmm <- as.factor(gmm_fit_final$classification)
#6. 查看最终各亚型人数分布
cat("=== 4 组亚型人数最终分布 ===\n")
print(table(nhanes_analytic_clustered$cluster_gmm, useNA = "ifany"))
