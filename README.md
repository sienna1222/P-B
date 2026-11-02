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
# Table1 人群基线表
#1. 加载必要的统计包
library(tableone)
library(dplyr)
#2. 定义你要展示在 Table 1 中的所有变量清单
#包含：人口学特征、生活方式、重新定义的既往疾病史、以及聚类核心/衰老指标
vars_to_compare <- c(
  "age", "gender", "race", "education", "poverty_ratio", 
  "bmi", "sbp", "dbp", "new_hyperten", "new_diabetes", 
  "phenoage0", "kdm0", "egfr"
)
#3. 确保分类变量（如性别、种族、教育、疾病史、聚类标签）在数据中被正确转换为因子类型
nhanes_analytic_clustered <- nhanes_analytic_clustered %>%
  mutate(
    gender = as.factor(gender),
    race = as.factor(race),
    education = as.factor(education),
    new_hyperten = as.factor(new_hyperten),
    new_diabetes = as.factor(new_diabetes),
    cluster_gmm = as.factor(cluster_gmm) # 确保聚类标签为因子
  )
#4. 明确指定哪些是分类变量（以便程序自动计算频数和构成比 %）
categorical_vars <- c("gender", "race", "education", "new_hyperten", "new_diabetes")
#5. 生成按 4 个聚类亚型（strata）分组的 Table 1 对象
table1_obj <- CreateTableOne(
  vars = vars_to_compare,
  strata = "cluster_gmm",
  data = nhanes_analytic_clustered,
  factorVars = categorical_vars,
  test = TRUE # 自动计算 p 值（连续变量用方差分析/秩和检验，分类变量用卡方检验）
)
#6. 在控制台中打印完整的 Table 1 结果（包含均值±标准差、频数/百分比及 p 值、SMD）
print(table1_obj, showAllLevels = TRUE, smd = TRUE)
#7. 导出为 CSV 文件，方便直接复制到 Word 或 Excel 中排版
table1_df <- print(table1_obj, printToggle = FALSE, showAllLevels = TRUE)
write.csv(table1_df, file = "Baseline_Table1_Final_Clusters.csv", row.names = TRUE)
cat("=== Table 1 已成功导出为 CSV 文件，可在当前工作目录下查看！===\n")
# K-M生存曲线
library(survminer)
km_plot_optimized <- ggsurvplot(
  km_fit,
  data = nhanes_analytic_clustered,
    #1. 核心显示优化
  pval = TRUE,             
  pval.coord = c(5, 0.82), # 将 p 值精确放置在图内部的合理坐标位置（可根据实际情况微调）
  conf.int = TRUE,         # 保留置信区间阴影（可设为 FALSE 让曲线更清爽）
    # 2. 纵轴范围缩放（核心优化：将纵轴截断显示为 0.75 到 1.00，放大曲线差异）
  ylim = c(0.75, 1.00),    
    # 3. 删失点与线条优化
  censv = TRUE,            # 是否显示删失点（若觉得太挤可设为 FALSE 不显示竖线）
  surv.median.line = "none", # 不画中位生存期虚线，保持画面干净
    # 4. 风险集表格（Risk Table）优化
  risk.table = TRUE,       
  risk.table.col = "strata",
  risk.table.height = 0.25,# 调整风险表的高度占比
  risk.table.y.text.col = TRUE, # 风险表左侧亚型名称带颜色
    # 5. 配色与标签美化
  palette = c("#2E9FDF", "#E7B800", "#FC4E07", "#2ca02c"), # 自定义经典 4 色
  xlab = "Follow-up Time (Months)", 
  ylab = "Survival Probability",
  legend.title = "Subtypes",
  legend.labs = c("Subtype 1", "Subtype 2", "Subtype 3", "Subtype 4"),
    # 6. 图例位置优化（放入图内右上角，节省版面）
  legend = c(0.85, 0.25), 
    ggtheme = theme_bw() + theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.background = element_rect(fill = alpha("white", 0.8))
  )
)
#打印优化后的图形
print(km_plot_optimized)
#导出为高分辨率 PDF（SCI 论文最爱用的矢量格式）
ggsave("KM_Survival_Curve_Optimized.pdf", plot = print(km_plot_optimized), width = 8, height = 7)
#导出为 300 DPI 的高清晰 TIFF 图片
tiff("KM_Survival_Curve_Optimized.tiff", width = 2400, height = 2100, res = 300)
print(km_plot_optimized)
dev.off()
