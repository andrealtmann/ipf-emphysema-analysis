library(lme4)
library(ggplot2)
library(gridExtra)

#analyze longitudinal datbase
filename <- "../data/Combined_RBHT_v2.csv"

fname.rd <- gsub(".csv",".RData", filename)

if (!file.exists(fname.rd)){
  message("running reshape_data.R")
  source("reshape_data.R")
} else {
  load(fname.rd)
}

years_since_ct <- long.df$days_since_ct / 365.25
long.df <- data.frame(long.df, years_since_ct)

## some univariate visualizations ##

#baseline FVC based on emph_groups#

long.ids <- unique(subset(long.lft, subset=max_days > 60)$ID)
dummy <- subset(sinfo, subset=is.element(Anonymised_CT_code, long.ids))

vars <- c("Age","Pk_yrs","Follow_uptodeath","FEV1","FEV1Pred","FVC","FVCPred","TLCOc_SB","dlco","CPI","abs_decline","rel_decline","TotalEmphave","ADMIXEDEMPHave","PUREEMPHave")
#dummy2 <- c()
#for(x in vars){
#  tmp <- data.frame(dummy$Anonymised_CT_code, dummy$emph_grp, x, dummy[,x])
#  colnames(tmp) <- c("ID","Emph","Variable","Value")
#  dummy2 <- rbind(dummy2, tmp)
#}
#ggplot(data = dummy2, aes(x=Emph,y=Value)) + geom_boxplot() + facet_wrap(~Variable,ncol=5)

plots <- list()
for(x in vars){
  fcol="white"
  pv <- anova(lm(as.formula(paste(x, "~ as.factor(emph_grp)")), data=dummy))[1,5]
  if (pv < 0.05){
    fcol="orange"
  }
  if (pv < 0.05 / length(vars)){
    fcol="red"
  }
  plots[[x]] <- ggplot(data = dummy, aes_string(x="emph_grp", y=x, group="emph_grp")) + geom_boxplot(fill=fcol)
}
grid.arrange(grobs=plots, ncol=4, nrow=4)

target <- paste(extract.lft, "x", sep=".")

lmer.formula <- paste(target, "- FVCPred", "~",  " as.factor(Scorers) + as.factor(Sex) + scaled_Age + Smokingneveris0 + as.factor(emph_grp) * years_since_ct + (1+years_since_ct|ID)")
lmer.formula2 <- paste(target, "- FVCPred", "~", " as.factor(Scorers) + as.factor(Sex) + scaled_Age + Smokingneveris0 + as.factor(emph_grp) + years_since_ct + (1+years_since_ct|ID)")

model <- lmer(as.formula(lmer.formula), data=long.df, subset=measures>1 & max_days > 90 )
model2 <- lmer(as.formula(lmer.formula2), data=long.df, subset=measures>1 & max_days > 90 )

anova(model, model2)
