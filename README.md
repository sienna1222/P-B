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
# 查看 4 个亚型在核心衰老/年龄指标上的均值，查找真正的“最优/最健康组”
nhanes_analytic_clustered %>%
  group_by(cluster_gmm) %>%
  summarise(
    N = n(),
    Mean_Age = mean(age, na.rm = TRUE),
    Mean_PhenoAge = mean(phenoage0, na.rm = TRUE),
    Mean_KDM = mean(kdm0, na.rm = TRUE),
    Death_Event_Rate = mean(status, na.rm = TRUE)
  )
  # 单因素Cox
library(survival)
#第一部分：数据预处理与单因素 Cox 回归
#1. 显式将 cluster_gmm 转换为因子，并将第 4 组（最健康/年轻组）设定为参照基准 (ref = "4", HR = 1.00)
nhanes_analytic_clustered$cluster_gmm <- relevel(factor(nhanes_analytic_clustered$cluster_gmm), ref = "4")
#2. 拟合单因素 Cox 比例风险回归模型
uni_cox_model <- coxph(
  Surv(permth_int, status) ~ cluster_gmm, 
  data = nhanes_analytic_clustered
)
#3. 提取详细结果并构建标准 SCI 论文格式的表格
uni_summary <- summary(uni_cox_model)
uni_hr <- exp(uni_summary$coefficients[, "coef"])
uni_ci <- exp(confint(uni_cox_model))
uni_pval <- uni_summary$coefficients[, "Pr(>|z|)"]
uni_results_table <- data.frame(
  Subtype = c("Subtype 1 vs Subtype 4", "Subtype 2 vs Subtype 4", "Subtype 3 vs Subtype 4"),
  HR = round(uni_hr, 3),
  Lower_95CI = round(uni_ci[, 1], 3),
  Upper_95CI = round(uni_ci[, 2], 3),
  P_value = signif(uni_pval, 3)
)
rownames(uni_results_table) <- NULL
cat("=== 1. 单因素 Cox 回归结果汇总表 ===\n")
print(uni_results_table)
#导出为表格文件
write.csv(uni_results_table, file = "Optimized_Univariate_Cox_Results.csv", row.names = FALSE)
#第二部分：等比例风险假定检验 (PH 假定检验)
#1. 运行严谨的统计学量化检验 (cox.zph)
ph_test <- cox.zph(uni_cox_model)
cat("\n=== 2. 等比例风险假定 (cox.zph) 统计检验结果 ===\n")
print(ph_test)
#2. 绘制发表级的 Schoenfeld 残差诊断图
#建议直接输出为高分辨率 PDF 矢量图供论文使用
pdf("Optimized_Schoenfeld_Residuals_Plot.pdf", width = 7, height = 6)
par(mfrow = c(1, 1), mar = c(4.5, 4.5, 2, 1))
plot(
  ph_test, 
  col = "black", 
  lwd = 1.5, 
  lty = c(1, 2, 2), # 实线为平滑拟合曲线，虚线为 95% 置信区间
  xlab = "Follow-up Time (Months)", 
  ylab = "Beta(t) for Cluster Subtypes"
)
abline(h = 0, col = "red", lty = 2, lwd = 1) # 添加 y=0 的参考基准线
title(main = "Schoenfeld Residuals Test for Proportional Hazards Assumption")
dev.off()
cat("\n=== 提示：Schoenfeld 诊断图已成功保存为矢量 PDF 文件！===\n")
# 绘制热图
library(dplyr)
library(tidyr)
#计算均值
#定义我们用于聚类和展示的核心变量
cluster_vars <- c("phenoage0", "kdm0", "egfr", "bmi", "sbp")
#计算各亚型在各核心指标上的均值
subtype_profile <- nhanes_analytic_clustered %>%
  group_by(cluster_gmm) %>%
  summarise(across(all_of(cluster_vars), ~ mean(.x, na.rm = TRUE)))
#对数据进行 Z-score 标准化（均值为0，方差为1）
subtype_profile_scaled <- subtype_profile %>%
  mutate(across(all_of(cluster_vars), ~ scale(.x)[, 1]))
print(subtype_profile_scaled)
#加载必要的精细化作图包
library(ggplot2)
library(dplyr)
library(tidyr)
#1. 确保数据准备就绪，并【强制设定因子的上下排列顺序】
#逻辑：将最健康的 Subtype 4 放在最下方，高危的 Subtype 2 放在最上方
subtype_profile_scaled$cluster_gmm <- factor(
  subtype_profile_scaled$cluster_gmm, 
  levels = c("2", "3", "1", "4") # ggplot2 默认纵轴自下而上，所以倒序排列可使 4 在最下、2 在最上
)
#2. 将宽数据转换为 ggplot2 专用的长数据格式 (Long format)
heatmap_long <- subtype_profile_scaled %>%
  pivot_longer(
    cols = all_of(c("phenoage0", "kdm0", "egfr", "bmi", "sbp")),
    names_to = "Variable",
    values_to = "Z_Score"
  )
#规范化变量显示名称
heatmap_long$Variable <- factor(
  heatmap_long$Variable,
  levels = c("phenoage0", "kdm0", "egfr", "bmi", "sbp"),
  labels = c("PhenoAge", "KDM Biological Age", "eGFR", "BMI", "Systolic BP")
)
heatmap_long$Subtype <- paste0("Subtype ", heatmap_long$cluster_gmm)
heatmap_long$Subtype <- factor(
  heatmap_long$Subtype, 
  levels = c("Subtype 2", "Subtype 3", "Subtype 1", "Subtype 4") # 确保视觉上由上至下符合高危->健康或健康->高危
)
#3. 使用 ggplot2 绘制高颜值、顶刊风格的精细化热图
p_heatmap_pro <- ggplot(heatmap_long, aes(x = Variable, y = Subtype, fill = Z_Score)) +
  #绘制带有精致白色网格边缘的方块
  geom_tile(color = "white", linewidth = 1.2, na.rm = FALSE) +
  #在方块中央精准标注标准化数值（保留两位小数）
  geom_text(aes(label = sprintf("%.2f", Z_Score)), color = "black", size = 4.5, fontface = "bold") +
  #采用高档的蓝-白-红发散配色 (RdBu 经典风格)
  scale_fill_gradient2(
    low = "#2166AC",      # 优质深蓝（代表低风险/负向偏离）
    mid = "#FFFFBF",      # 柔和过渡白黄色
    high = "#B2182B",     # 优质深红（代表高风险/正向偏离）
    midpoint = 0,
    limit = c(-1.5, 1.5), 
    name = "Z-Score"
  ) +
  #极简且专业的学术主题排版
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, face = "bold", color = "black", size = 12),
    axis.text.y = element_text(face = "bold", color = "black", size = 12),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16, margin = margin(b = 15)),
    legend.title = element_text(face = "bold", size = 11),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(
    title = "Clinical and Aging Profiles Across GMM Subtypes"
  )
#4. 在控制台打印预览
print(p_heatmap_pro)
#5. 导出为高分辨率出版级 PDF 矢量图
ggsave("Publication_Quality_Heatmap_Reordered.pdf", plot = p_heatmap_pro, width = 8, height = 5, dpi = 300)
cat("\n=== 提示：纵轴排序修正后的顶刊热图已成功导出为 PDF 文件！===\n")
