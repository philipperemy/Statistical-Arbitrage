#http://quanttrader.info/public/testForCoint.html

library(zoo)
gld <- read.csv("http://ichart.finance.yahoo.com/table.csv?s=BNP.PA&ignore=.csv", stringsAsFactors=F)
gdx <- read.csv("http://ichart.finance.yahoo.com/table.csv?s=GLE.PA&ignore=.csv", stringsAsFactors=F)

gld_dates <- as.Date(gld[,1])
gdx_dates <- as.Date(gdx[,1])

gld <- zoo(gld[,7], gld_dates)
gdx <- zoo(gdx[,7], gdx_dates)

t.zoo <- merge(gld, gdx, all=FALSE)

t <- as.data.frame(t.zoo)

cat("Date range is", format(start(t.zoo)), "to", format(end(t.zoo)), "\n")

m <- lm(gld ~ gdx + 0, data=t)
beta <- coef(m)[1]
cat("Assumed hedge ratio is", beta, "\n")

sprd <- t$gld - beta*t$gdx

library(tseries)

sprd2 <- sprd[100:length(sprd)]
ht <- adf.test(sprd2, alternative="stationary", k=0)
cat("ADF p-value is", ht$p.value, "\n")

if (ht$p.value < 0.05) {
  cat("The spread is likely mean-reverting.\n")
} else {
  cat("The spread is not mean-reverting.\n")
}