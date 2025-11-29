
library(ggforce)
library(ggplot2)

df1 = readxl::read_xlsx("1_BMI_K=4_migraine_direct.xlsx")

df1$exposure_outcome = factor(df1$exposure_outcome, levels =  c("1-year-old BMI -> 8-year-old BMI", "1-year-old BMI -> Adult BMI", "8-year-old BMI -> Adult BMI",  "1-year-old BMI -> Migraine, total",  
                                                                "8-year-old BMI -> Migraine, total", "Adult BMI -> Migraine, total" )   )

df1$color_group <- ifelse(df1$OR_lower > 1 & df1$OR_upper > 1, "red",
                          ifelse(df1$OR_lower < 1 & df1$OR_upper < 1, "blue", "black"))


p1 = ggplot(df1, aes(y = exposure_outcome, x = OR)) +
  geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper, color = color_group), height = 0.2, size = 0.6) +
  geom_point(aes(x = OR, color = color_group), size = 2) +  # Add points for mean
  scale_color_manual(values = c("black" = "black", "red" = "#FF6B6B", "blue" = "#6A9BFF")) +
  guides(color = "none")+
  theme_bw() +
  theme(text = element_text(family = "sans"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y.left = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(angle = 0),
        plot.title = element_text(size = 14, face = "bold"))+
  geom_vline(xintercept = 1, size = 0.7)+
  ggtitle("Direct Effects (OR)")
p1

ggsave(
  filename = "1_BMI_K=4_migraine_direct_plot.tiff", 
  plot = p1,                            
  device = "tiff",                      
  dpi = 600,                            
  width = 5,                            
  height = 3,                           
  units = "in"                          
)


