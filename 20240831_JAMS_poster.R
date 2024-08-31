#------------------------------------------------------------
#This R script was presented as a poster at the 77th Annual Meeting of the Japanese Association for Mathematical Sociology on August 31, 2024. 
#The title was "Causal Effects of Enrollment Shortfalls on Financial Performance in Private Universities: An Application of G-Methods."
#If you need a translation of the following description, which is in Japanese, please feel free to contact me, and I will provide you with an English translation.
#------------------------------------------------------------

# ------------------------------------------------------------
# 準備
# ------------------------------------------------------------
# Rのバージョン確認
R.version

# ワークスペースのクリア: オブジェクト一括削除
rm(list = ls())

# データ取り込み
df <- read.csv(file.choose())
head(df, 5)

#------------------------------------------------------------
#パッケージの読み込み
#------------------------------------------------------------
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)
if (!require("mgcv")) install.packages("mgcv")
library(mgcv)
if (!require("viridis")) install.packages("viridis")
library(viridis)
if (!require("patchwork")) install.packages("patchwork")
library(patchwork)
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)
if (!require("tidyr")) install.packages("tidyr")
library(tidyr)
if (!require("RColorBrewer")) install.packages("RColorBrewer")
library(RColorBrewer)

# ------------------------------------------------------------
# データの前処理
# ------------------------------------------------------------
# fillrの計算
df <- df %>%
  mutate(fillr = enrolment / capa * 100)  # fillrを百分率で計算

# 時間依存性変数の作成
df <- df %>%
  group_by(id) %>%
  arrange(wave) %>%
  mutate(
    lag_fillr = lag(fillr),
    lag_outcome = lag(outcome),
    lag_students = lag(students),
    lag_status = lag(status),
    lag_confounder = lag(confounder)
  ) %>%
  ungroup()

# 完全データのみを選択
df_complete <- df %>%
  select(id, wave, fillr, lag_fillr, lag_outcome, lag_students, lag_status, lag_confounder,
         scienced, hospitald, womend, cityd, foundation, outcome, confounder) %>%
  drop_na()

# 事業活動収支差額の対数変換（負の値に対応するため，対称的対数変換（Symmetric Log Transform）を実施）
df_complete$log_outcome <- sign(df_complete$outcome) * log(abs(df_complete$outcome) + 1)

# ------------------------------------------------------------
# 散布図の作成
# ------------------------------------------------------------
p_base <- ggplot(df_complete, aes(x = fillr, y = outcome)) +
  geom_point(alpha = 0.3) +  # 透明度を設定して重なりを表現
  geom_smooth(method = "loess", color = "red") +  # 非線形の傾向線を追加
  labs(x = "入学定員充足率 (%)", 
       y = "事業活動収支差額",
       title = "入学定員充足率と経常収支差額の関係") +
  theme_minimal()

# 対数変換したoutcomeでも散布図を作成
p_log_base <- ggplot(df_complete, aes(x = fillr, y = log_outcome)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "red") +
  labs(x = "入学定員充足率 (%)", 
       y = "事業活動収支差額（対数変換）",
       title = "入学定員充足率と事業活動収支差額（対数変換）の関係") +
  theme_minimal()

# プロットの表示
print(p_base)
print(p_log_base)

# ------------------------------------------------------------
# IPWの計算(Generalized Additive Modelを使用)
# ------------------------------------------------------------
# waveを数値として扱う（2013-2020の8年分）
df_complete$wave_numeric <- as.numeric(as.character(df_complete$wave))

ps_model <- gam(
  fillr ~ s(lag_fillr) + s(lag_outcome) + s(lag_students) + s(lag_status) + s(lag_confounder) +
    scienced + hospitald + womend + cityd + foundation + 
    s(wave_numeric, bs = "cr", k = 8),  # 8年分のデータなので、k=8を使用
  family = gaussian(),
  data = df_complete
)

# IPWの計算
df_complete$ps <- predict(ps_model, type = "response")

# 分母が0になるのを防ぐため、小さな値を加える
epsilon <- 1e-8
df_complete$ipw <- 1 / (df_complete$ps + epsilon)

# 安定化重みの計算
stabilized_model <- gam(
  fillr ~ s(wave_numeric, bs = "cr", k = 8),
  family = gaussian(),
  data = df_complete
)
df_complete$ps_stabilized <- predict(stabilized_model, type = "response")
df_complete$stabilized_ipw <- df_complete$ipw * (df_complete$ps_stabilized + epsilon)

# 重みの確認
summary(df_complete$stabilized_ipw)
hist(df_complete$stabilized_ipw, main = "Distribution of Stabilized IPW", xlab = "Weight")

# 重みの切り詰め
lower_bound <- quantile(df_complete$stabilized_ipw, 0.01)
upper_bound <- quantile(df_complete$stabilized_ipw, 0.99)
df_complete$trimmed_ipw <- pmax(pmin(df_complete$stabilized_ipw, upper_bound), lower_bound)

# 切り詰め後の重みの確認
summary(df_complete$trimmed_ipw)
hist(df_complete$trimmed_ipw, main = "Distribution of Trimmed IPW", xlab = "Weight")

# ------------------------------------------------------------
# MSMの適用(Generalized Additive Modelを使用)
# ------------------------------------------------------------
# GAMを使用したMSMの適用
msm_gam_log <- gam(
  log_outcome ~ s(fillr, bs = "cr"),
  weights = trimmed_ipw,
  data = df_complete
)

summary(msm_gam_log)

# ------------------------------------------------------------
# MSMの適用(Generalized Additive Modelを使用)：結果の可視化
# ------------------------------------------------------------
# GAMモデルから予測データを取得
pred_data <- predict(msm_gam_log, type = "terms", se = TRUE)

# データフレームを作成
plot_data <- data.frame(
  fillr = df_complete$fillr,
  effect = pred_data$fit[, "s(fillr)"],
  se = pred_data$se.fit[, "s(fillr)"]
)

# ggplotでグラフを作成
ggplot(plot_data, aes(x = fillr, y = effect)) +
  geom_ribbon(aes(ymin = effect - 1.96 * se, ymax = effect + 1.96 * se), 
              alpha = 0.2) +
  geom_line() +
  geom_rug(sides = "b", alpha = 0.1) +
  labs(
    x = "入学定員充足率（%）",
    y = "経常収支差額への効果（対数スケール）",
    title = "入学定員充足率が経常収支差額に与える非線形効果",
    subtitle = "信頼区間（95%）付き"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  scale_x_continuous(breaks = seq(0, max(df_complete$fillr), by = 50)) +
  coord_cartesian(ylim = c(min(plot_data$effect - 1.96 * plot_data$se),
                           max(plot_data$effect + 1.96 * plot_data$se)))

# fillrが120%以下のデータのみを使用
plot_data_filtered <- plot_data[plot_data$fillr <= 120, ]

# ggplotでグラフを作成
ggplot(plot_data_filtered, aes(x = fillr, y = effect)) +
  geom_ribbon(aes(ymin = effect - 1.96 * se, ymax = effect + 1.96 * se), 
              alpha = 0.2) +
  geom_line() +
  geom_rug(sides = "b", alpha = 0.1) +
  labs(
    x = "入学定員充足率（%）",
    y = "事業活動収支差額への効果（対数スケール）"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  scale_x_continuous(breaks = seq(0, 120, by = 10), limits = c(0, 120)) +
  coord_cartesian(ylim = c(min(plot_data_filtered$effect - 1.96 * plot_data_filtered$se),
                           max(plot_data_filtered$effect + 1.96 * plot_data_filtered$se)))

# ------------------------------------------------------------
# 交互作用を含むMSMの適用（GAM,cubic regression splineを使用）
# ------------------------------------------------------------
# 交互作用を含むGAMモデルの適用（対数変換後）
msm_gam_interaction_log <- gam(
  log_outcome ~ te(fillr, lag_confounder, bs = "cr"),
  weights = trimmed_ipw,
  data = df_complete
)

# モデルのサマリーを表示
summary(msm_gam_interaction_log)

# ------------------------------------------------------------
# 交互作用を含むMSMの適用（GAM,cubic regression splineを使用）：結果の可視化
# ------------------------------------------------------------
# データの実際の範囲を確認
min_fillr <- min(df_complete$fillr)
max_fillr <- 120  # 上限を120%に設定

# 予測データの生成（下限は実際のデータ範囲、上限は120%に制限）
new_data <- expand.grid(
  fillr = seq(min_fillr, max_fillr, length.out = 100),
  lag_confounder = c(-1, 0, 1)
)

# モデルを使用して予測
predictions <- predict(msm_gam_interaction_log, newdata = new_data, se.fit = TRUE)
new_data$predicted <- predictions$fit
new_data$se <- predictions$se.fit

# 信頼区間の計算
new_data$lower <- new_data$predicted - 1.96 * new_data$se
new_data$upper <- new_data$predicted + 1.96 * new_data$se

# lag_confounderにラベルを付ける
new_data$lag_confounder_label <- factor(new_data$lag_confounder,
                                        levels = c(-1, 0, 1),
                                        labels = c("減少", "変化なし", "増加"))

# RdYlBuパレットから色を選択
colors <- brewer.pal(3, "RdYlBu")

# プロットの作成（範囲を制限）
p_log <- ggplot(new_data, aes(x = fillr, y = predicted, color = lag_confounder_label)) +
  geom_line(aes(linetype = lag_confounder_label), linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = lag_confounder_label), alpha = 0.1) +
  labs(x = "入学定員充足率 (%)", 
       y = "事業活動収支差額への効果（対数スケール）",
       color = "前年度からの\n入学定員数変化",
       fill = "前年度からの\n入学定員数変化",
       linetype = "前年度からの\n入学定員数変化") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_x_continuous(limits = c(min_fillr, max_fillr), 
                     breaks = seq(ceiling(min_fillr/10)*10, 120, by = 10)) +
  theme_minimal() +
  theme(legend.position = "bottom") 

print(p_log)