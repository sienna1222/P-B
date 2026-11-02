# 数据载入
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
library(dplyr)
library(mclust)
# 数据排纳
#步骤 1：年龄筛选 (纳入 >= 45 岁中老年人群)
step1_data <- nhanes_work %>% filter(!is.na(age) & age >= 45)
N1 <- nrow(step1_data)
cat("步骤 1 (筛选年龄 >= 45 岁): 剔除", n_current - N1, "人 | 剩余:", N1, "人\n") # 已修复 n_step1 为 N1
n_current <- N1
#步骤 2：排除孕妇 (pregnant == 1)
step2_data <- step1_data %>% filter(pregnant != 1 | is.na(pregnant))
N2 <- nrow(step2_data)
cat("步骤 2 (排除孕妇): 剔除", n_current - N2, "人 | 剩余:", N2, "人\n")
n_current <- N2
#步骤 3：剔除核心原生变量、随访数据及聚类变量的缺失者
#把 cluster_vars（phenoage0, kdm0, egfr, bmi, sbp）全部纳入过滤，确保零 NA
step3_data <- step2_data %>% filter(
  !is.na(gender),
  !is.na(creat),
  !is.na(egfr),
  !is.na(phenoage0),
  !is.na(kdm0),
  !is.na(status),
  !is.na(permth_int),
  !is.na(bmi),
  !is.na(sbp)
)
N3 <- nrow(step3_data)
cat("步骤 3 (剔除核心变量及随访/聚类指标缺失者): 剔除", n_current - N3, "人 | 剩余:", N3, "人\n")
n_current <- N3
#步骤 4：剔除极端离群值 (eGFR < 15)
nhanes_analytic <- step3_data %>% filter(egfr >= 15)
N4 <- nrow(nhanes_analytic)
cat("步骤 4 (剔除 eGFR < 15 极端值): 剔除", n_current - N4, "人 | 最终分析样本:", N4, "人\n")
cat("========================\n")

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
