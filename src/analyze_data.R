library(lme4)
library(ggplot2)
library(gridExtra)
library(survival)

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

#raw score
#target <- paste(extract.lft, "x", sep=".")
#abs difference to t0
target <- paste(extract.lft, ".x - FVC.y", sep="")
#rel difference to t0
#target <- paste("(", extract.lft, ".x - FVC.y",") / FVC.y", sep="")


#### univariate feature screen ####
svars <- colnames(sinfo)[11:40]
svars <- setdiff(svars, c("CTdate","scaled_Age"))
##renamve FVC
svars[svars == "FVC"] <- "FVC.y"
#svars <- gsub("emph_grp", "as.factor(emph_grp)", svars)

#with lme
res.lme <- t(sapply(svars, function(vvv){
  message(vvv)
  base.formula <- paste(target, "~",  "FVC.y + as.factor(Scorers) + as.factor(Sex) + scaled_Age + Smokingneveris0 + Pk_yrs + years_since_ct + (1+years_since_ct|ID)")
  uni.formula  <- paste(target, "~",  "FVC.y + as.factor(Scorers) + as.factor(Sex) + scaled_Age + Smokingneveris0 + Pk_yrs + years_since_ct", "+",  vvv, "+ (1+years_since_ct|ID)")
  int.formula  <- paste(target, "~",  "FVC.y + as.factor(Scorers) + as.factor(Sex) + scaled_Age + Smokingneveris0 + Pk_yrs + years_since_ct", "*",  vvv, "+ (1+years_since_ct|ID)")

  use.df <- subset(long.df, subset=measures>1 & max_days > 60)
  idx <- !is.na(use.df[,vvv])

  base.model   <- lmer(as.formula(base.formula), data=use.df, subset= idx)
  uni.model   <- lmer(as.formula(uni.formula), data=use.df, subset=idx)
  int.model   <- lmer(as.formula(int.formula), data=use.df, subset=idx)

  p1 <- anova(base.model, uni.model)[2,8]
  p2 <- anova(uni.model, int.model)[2,8]
  p3 <- anova(base.model, int.model)[2,8]
  return(c(p1,p2,p3))
}))

#with survival
svars2 <- gsub("FVC.y","FVC",svars)
res.cox <- t(sapply(svars2, function(vvv){
  message(vvv)
  #use.df <- subset(long.df, subset=measures>1 & max_days > 60)
  use.df <- sinfo
  idx <- !is.na(use.df[,vvv])
  mysurv <- Surv(use.df$Follow_uptodeath, use.df$Dead)

  base.formula <- paste("mysurv ~ FVC + as.factor(Scorers) + as.factor(Sex) + scaled_Age + Smokingneveris0 + Pk_yrs")
  uni.formula  <- paste("mysurv ~ FVC + as.factor(Scorers) + as.factor(Sex) + scaled_Age + Smokingneveris0 + Pk_yrs + ", vvv)

  base.cox <- coxph(as.formula(base.formula), data=use.df, subset=idx)
  uni.cox <- coxph(as.formula(uni.formula), data=use.df, subset=idx)


  pv <- anova(base.cox, uni.cox)[2,4]
  zs <- summary(uni.cox)$coeff[vvv,4]

  return( c(pv, zs))
}))

### addressing the main question ###

lmer.formula <- paste(target, "- FVC.y", "~",  " as.factor(Scorers) + as.factor(Sex) + scaled_Age + Smokingneveris0 + as.factor(emph_grp) * years_since_ct + (1+years_since_ct|ID)")
lmer.formula2 <- paste(target, "- FVC.y", "~", " as.factor(Scorers) + as.factor(Sex) + scaled_Age + Smokingneveris0 + as.factor(emph_grp) + years_since_ct + (1+years_since_ct|ID)")
model <- lmer(as.formula(lmer.formula), data=long.df, subset=measures>1 & max_days > 60)
model2 <- lmer(as.formula(lmer.formula2), data=long.df, subset=measures>1 & max_days > 60)
anova(model, model2)

#with cox regression
emphvar <- c("PUREEMPHave","ADMIXEDEMPHave","TotalEmphave","Emph_presence","emph_grp","Emphover5","Emphover10","Emphover15")
use.df <- subset(sinfo, subset=dlco>0)
mysurv <- Surv(use.df$Follow_uptodeath, use.df$Dead)
base.formula <- paste("mysurv ~ FVC + GAPindex + as.factor(Scorers) + as.factor(Sex) + scaled_Age + Smokingneveris0 + Pk_yrs")
sapply(emphvar, function(ev){
  uni.formula  <- paste(base.formula, "+", ev)
  idx <- !is.na(use.df[,ev])
  base.cox <- coxph(as.formula(base.formula), data=use.df, subset=idx)
  uni.cox <- coxph(as.formula(uni.formula), data=use.df, subset=idx)

  pv <- anova(base.cox, uni.cox)[2,4]
  return(pv)
})
