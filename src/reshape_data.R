
filename <- "../data/Combined_RBHT_v2.csv"

mydata <- read.csv(filename)

#contains also most lung function values from time-point one
sinfo <- mydata[,1:38]

emph_grp <- apply(sinfo,1,function(x){

  res <- 0
  if (x["Emph_presence"]==1){
    res <- 1
  }
  if (x["Emphover5"]==1){
    res <- 2
  }
  if (x["Emphover10"]==1){
    res <- 3
  }
  if (x["Emphover15"]==1){
    res <- 4
  }
  return(res)
})

sinfo <- data.frame(sinfo, scale(sinfo[,"Age"]), emph_grp)
colnames(sinfo)[39] <- "scaled_Age"


#get longitudinal data
lftdates <- paste("LFT", paste(1:18), sep="")
ctdate   <- "CTdate"
extract.lft <- "FVC"

long.lft <- c()
dev.null <- apply(mydata,1,function(x){
  my.ct.date <- as.Date(paste(x[ctdate]), format="%m/%d/%y")
  res <- c()
  for (i in 1:18){
    d1 <- lftdates[i]
    d2 <- paste(extract.lft, i, sep="")
    if (x[d1] != "" & !is.na(x[d2])){
      d1.date <- as.Date(paste(x[d1]), format="%m/%d/%y")
      res <- rbind(res, c(as.integer(x[1]), d1.date-my.ct.date, as.double(x[d2])))
    }
  }
  nr <- nrow(res)
  md <- res[nr,2]
  res <- cbind(res, nr, md)
  long.lft <<- rbind(long.lft, res)
  return(res)
})
colnames(long.lft) <- c("ID","days_since_ct", extract.lft, "measures", "max_days")
long.lft <- data.frame(long.lft)

##derive annualized decline##
lm.formula <- paste(extract.lft, "~ days_since_ct")
an.decline <- t(sapply(unique(long.lft$ID), function(id){
  message(id)
  tmp <- subset(long.lft, subset=ID==id)
  if (nrow(tmp)==1){
    return(c(id, rep(NA,2)))
  }
  mmm <- lm(as.formula(lm.formula), data=tmp)
  abs.decline <- mmm$coeff[2] * 365.25
  rel.decline <- abs.decline / mmm$coeff[1]
  return(c(id, abs.decline, rel.decline))
}))
colnames(an.decline) <- c("ID", "abs_decline", "rel_decline")

sinfo <- merge(sinfo, an.decline, by.x="Anonymised_CT_code", by.y="ID", all.x=T)

long.df <- merge(long.lft, sinfo, by.x="ID", by.y="Anonymised_CT_code")

save(sinfo, long.lft, long.df, extract.lft, file=gsub(".csv",".RData",filename))
