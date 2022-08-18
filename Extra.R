getwd()
setwd("/Users/chensisi/Documents/RNAseq/4_Cell_deconvolution/")

load("celldecon.RData")

LGG.merge$group <- NA
High <- which(LGG.merge$T.cell.CD4. > quantile(LGG.merge$T.cell.CD4., 0.75))
Low<- which(LGG.merge$ENST00000620340.4 < quantile(LGG.merge$ENST00000620340.4, 0.25))
overlap <- intersect(High, Low)
High <- High[-which(High %in% overlap)]
Low <- Low[-which(Low %in% overlap)]
LGG.merge$group[High] <- "high"
LGG.merge$group[Low] <- "low"
LGG.merge$group[which(LGG.merge$ENST00000620340.4 < -8)] <- NA
LGG.merge$group[which(LGG.merge$T.cell.CD4. > 0.1 & LGG.merge$ENST00000620340.4 > -0.5)] <- "high"

table(LGG.merge$group)

p <- ggplot(LGG.merge, aes(ENST00000620340.4, T.cell.CD4.)) + 
  geom_point() + stat_smooth(method = "nls", formula = 'y~a*exp(b*x)',method.args = list(start=c(a=0.1646, b=9.5e-8)), se=FALSE) + 
  stat_cor() + theme_bw() + my_theme + 
  labs(x = "Isoform 1 expression", y = "Proportion of CD4+ T cells")
q <- p + geom_point(aes(ENST00000620340.4, T.cell.CD4., color = group)) 
plot(q)
# ggsave(q, file = "fig4d.pdf",
#       height = 15, width = 15, units = "cm")

test <- LGG.merge %>% filter(group == "high"| group == "low")
LGG.os <- Surv(test$OS.time, test$OS)

# Survival analysis
summary(coxph(LGG.os ~ test$T.cell.CD4.))
test <- test %>% mutate(CD4.high = ifelse(group == "high", 1, 0))
CD4.high <- as.numeric(test$CD4.high)

fit <- survfit(LGG.os ~ test$CD4.high)
legend.names <- c(paste0("CD4 T cell","-low"), paste0("CD4 T cell","-high")) # IF / ratio 
ggsurvplot(fit, data= LGG.os,
           pval = TRUE, pval.method = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           risk.table.col = "strata", # Change risk table color by group 
           #surv.median.line = "hv", # Specify median survival
           legend.title = "",
           legend.labs= legend.names,
           palette = c("#00BFC4","#F8766D"),
           xlab= "OS time (days)",
           ylab= "Proportion survival rate",
           ggtheme = theme_bw() + my_theme)








