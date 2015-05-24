#http://gekkoquant.com/2012/10/21/statistical-arbitrage-correlation-vs-cointegration/

ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}

CCE <- read.csv("CCE.csv")
PEP <- read.csv("PEP.csv")

C <- rev(CCE[,5])
P <- rev(PEP[,5])

cor(C, P)

plot(C, type="l")
plot(P, type="l")

lm.1 <- lm(C ~ P)
summary(lm.1)
beta <- 0.518715
alpha <- -7.156312

S <- C - beta*P - alpha
plot(S, type="l")

plot(C, type="l", col="Blue")
lines(beta*P + alpha, type="l", col="Red")

# Middle Band = 20-day simple moving average (SMA)
# Upper Band = 20-day SMA + (20-day standard deviation of price x 2) 
# Lower Band = 20-day SMA - (20-day standard deviation of price x 2)
plot(S, type="l")
MA_S <- ma(S, 50)
lines(MA_S, type="l", col="Red")
lines(MA_S + 2*sd, type="l", col="Blue") #false for SD calculation
lines(MA_S - 2*sd, type="l", col="Blue")
mean(S)
sd(S)

library(tseries)
adf.test(S,k=0)

#not stationary

#I(1)
diff_C <- diff(C)
diff_P <- diff(P)

adf.test(diff_C,k=0)
adf.test(diff_P,k=0)

library(moments)


#### FX ####

plot(EUR_GBP_Week1$RateAsk, type="l")

EURGBP_ASK <- EUR_GBP_Week1$RateAsk
EURUSD_ASK <- EUR_USD_Week1$RateAsk

EURUSD_ASK = EURUSD_ASK[1:343636]

lm.3 <- lm(EURGBP_ASK ~ EURUSD_ASK)
summary(lm.3)

Spr <- EURGBP_ASK - 0.002452 * EURUSD_ASK - 0.748242
plot(Spr, type="l")
adf.test(Spr,k=0)
