library("quantmod")
library("PerformanceAnalytics")
library("fUnitRoots")
library("tseries")

#Specify dates for downloading data, training models and running simulation
hedgeTrainingStartDate = as.Date("2008-01-01") #Start date for training the hedge ratio
hedgeTrainingEndDate = as.Date("2010-01-01") #End date for training the hedge ratio

symbolLst<-c("RDS-A","RDS-B")
title<-c("Royal Dutch Shell A vs B Shares")


### SECTION 1 - Download Data & Calculate Returns ###
#Download the data
symbolData <- new.env() #Make a new environment for quantmod to store data in
getSymbols(symbolLst, env = symbolData, src = "yahoo", from = hedgeTrainingStartDate, to=hedgeTrainingEndDate)

stockPair <- list(
  a = coredata(Cl(eval(parse(text=paste("symbolData$\"",symbolLst[1],"\"",sep="")))))   #Stock A
  ,b = coredata(Cl(eval(parse(text=paste("symbolData$\"",symbolLst[2],"\"",sep=""))))) #Stock B
  ,name=title)




testForCointegration <- function(stockPairs){
  #Pass in a pair of stocks and do the necessary checks to see if it is cointegrated
  
  #Plot the pairs
  
  dev.new()
  par(mfrow=c(2,2))
  print(c((0.99*min(rbind(stockPairs$a,stockPairs$b))),(1.01*max(rbind(stockPairs$a,stockPairs$b)))))
  plot(stockPairs$a, main=stockPairs$name, ylab="Price", type="l", col="red",ylim=c((0.99*min(rbind(stockPairs$a,stockPairs$b))),(1.01*max(rbind(stockPairs$a,stockPairs$b)))))
  lines(stockPairs$b, col="blue")
  print(length(stockPairs$a))
  print(length(stockPairs$b))
  #Step 1: Calculate the daily returns
  dailyRet.a <- na.omit((Delt(stockPairs$a,type="arithmetic")))
  dailyRet.b <- na.omit((Delt(stockPairs$b,type="arithmetic")))
  dailyRet.a <- dailyRet.a[is.finite(dailyRet.a)] #Strip out any Infs (first ret is Inf)
  dailyRet.b <- dailyRet.b[is.finite(dailyRet.b)]
  
  print(length(dailyRet.a))
  print(length(dailyRet.b))
  #Step 2: Regress the daily returns onto each other
  #Regression finds BETA and C in the linear regression retA = BETA * retB + C
  regression <- lm(dailyRet.a ~ dailyRet.b + 0)
  beta <- coef(regression)[1]
  print(paste("The beta or Hedge Ratio is: ",beta,sep=""))
  plot(x=dailyRet.b,y=dailyRet.a,type="p",main="Regression of RETURNS for Stock A & B") #Plot the daily returns
  lines(x=dailyRet.b,y=(dailyRet.b*beta),col="blue")#Plot in linear line we used in the regression
  
  
  #Step 3: Use the regression co-efficients to generate the spread
  spread <- stockPairs$a - beta*stockPairs$b #Could actually just use the residual form the regression its the same thing
  spreadRet <- Delt(spread,type="arithmetic")
  spreadRet <- na.omit(spreadRet)
  #spreadRet[!is.na(spreadRet)]
  plot((spreadRet), type="l",main="Spread Returns") #Plot the cumulative sum of the spread
  plot(spread, type="l",main="Spread Actual") #Plot the cumulative sum of the spread
  #For a cointegrated spread the cumsum should not deviate very far from 0
  #For a none-cointegrated spread the cumsum will likely show some trending characteristics
  
  #Step 4: Use the ADF to test if the spread is stationary
  #can use tSeries library
  adfResults <- adf.test((spread),k=0,alternative="stationary")
  
  print(adfResults)
  if(adfResults$p.value <= 0.05){
    print(paste("The spread is likely Cointegrated with a pvalue of ",adfResults$p.value,sep=""))
  } else {
    print(paste("The spread is likely NOT Cointegrated with a pvalue of ",adfResults$p.value,sep=""))
  }
  
}


testForCointegration(stockPair)
