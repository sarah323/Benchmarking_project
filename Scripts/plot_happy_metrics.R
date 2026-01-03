# Simple plotting script for hap.py metrics on HG002 chr22

# Load required package
library(ggplot2)

# Helper function to read one summary file and keep PASS rows only
read_happy_summary <- function(file, depth_label) {
  df <- read.csv(file, header = TRUE)
  df <- subset(df, Filter == "PASS" & Type %in% c("SNP", "INDEL"))
  df$Depth <- depth_label
  return(df)
}

# Read the four coverage levels
d2  <- read_happy_summary("HG002.chr22_2x.summary.csv",  "2x")
d10 <- read_happy_summary("HG002.chr22_10x.summary.csv", "10x")
d40 <- read_happy_summary("HG002.chr22_40x.summary.csv", "40x")
d80 <- read_happy_summary("HG002.chr22_80x.summary.csv", "80x")

all_df <- rbind(d2, d10, d40, d80)

# Make Depth an ordered factor so it appears in the right order on the x-axis
all_df$Depth <- factor(all_df$Depth, levels = c("2x", "10x", "40x", "80x"))

#### 1) F1-score vs depth (SNPs vs INDELs) ####
p_f1 <- ggplot(all_df, aes(x = Depth, y = METRIC.F1_Score,
                           group = Type, linetype = Type)) +
  geom_line() +
  geom_point() +
  labs(title = "F1-score across coverage depth (chr22 HG002)",
       y = "F1-score", x = "Coverage depth") +
  theme_bw()

#### 2) Recall vs depth ####
p_recall <- ggplot(all_df, aes(x = Depth, y = METRIC.Recall,
                               group = Type, linetype = Type)) +
  geom_line() +
  geom_point() +
  labs(title = "Recall across coverage depth (chr22 HG002)",
       y = "Recall (Sensitivity)", x = "Coverage depth") +
  theme_bw()

#### 3) Precision vs depth ####
p_prec <- ggplot(all_df, aes(x = Depth, y = METRIC.Precision,
                             group = Type, linetype = Type)) +
  geom_line() +
  geom_point() +
  labs(title = "Precision across coverage depth (chr22 HG002)",
       y = "Precision (PPV)", x = "Coverage depth") +
  theme_bw()

# Save each plot as a PNG file (for the report)
ggsave("hap_py_F1_depth.png",      p_f1,     width = 6, height = 4, dpi = 300)
ggsave("hap_py_Recall_depth.png",  p_recall, width = 6, height = 4, dpi = 300)
ggsave("hap_py_Precision_depth.png", p_prec, width = 6, height = 4, dpi = 300)


