library(dplyr)
# 1. 载入原始数据集
data(NHANES4)
# 2. 定义 CKD-EPI 2021（无种族修正版）计算函数
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
# 3. 执行排纳流程
n_total <- nrow(NHANES4)
nhanes_analytic <- NHANES4 %>%
  # 新增 egfr 变量
  mutate(egfr = calc_egfr_2021(scr = creat, age = age, sex = gender)) %>%
  # 步骤 1：年龄筛选 (纳入 >= 45 岁成年人)
  filter(!is.na(age) & age >= 45) %>%
  { n_step1 <- nrow(.); message("Step 1 (Age >= 45): ", n_step1); . } %>%
  # 步骤 2：排除孕妇 (pregnant == 1)
  filter(pregnant != 1 | is.na(pregnant)) %>%
  { n_step2 <- nrow(.); message("Step 2 (Exclude Pregnant): ", n_step2); . } %>%
  # 步骤 3：剔除核心原生变量缺失者
  filter(
    !is.na(gender),
    !is.na(creat),
    !is.na(egfr),
    !is.na(phenoage0),
    !is.na(kdm0)
  ) %>%
  { n_step3 <- nrow(.); message("Step 3 (Complete Core Data): ", n_step3); . } %>%
  # 步骤 4：剔除极端离群值 (eGFR < 15)
  filter(egfr >= 15) %>%
  { n_step4 <- nrow(.); message("Step 4 (Final Analytic Sample): ", n_step4); . }
# 查看最终分析样本量
cat("原始总样本量:", n_total, "\n")
cat("最终纳入分析的样本量:", nrow(nhanes_analytic), "\n")
