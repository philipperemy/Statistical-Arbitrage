rho = 0.91;
sigma = 1;
beta = 0.5; 
steps = 10000;
nu = 3;

x = numeric(steps);
y = numeric(steps);
x[1] = rnorm(1) * sqrt(sigma^2/(1-rho^2));      

for (i in 1:steps) {
  x[i+1] = rho  * x[i]  + sigma*rnorm(1);
  y[i]   = beta * exp(0.5*x[i]) *rt(1, nu); #rnorm(1)
}


plot(y, type="l")

acf(y)
acf(y^2)

library(MASS)
truehist(y)
qqnorm(y, main = " ")
qqline(y)
hist(y)

install.packages("ghyp")
library(ghyp)


fit <- stepAIC.ghyp(y, dist = c("ghyp", "hyp", "NIG", "VG", "t", "gauss"),symmetric = NULL, control = list(maxit = 5000),na.rm = T, silent = TRUE)

#fit <- stepAIC.ghyp(y, dist = c("t"),symmetric = NULL, control = list(maxit = 5000),na.rm = T, silent = TRUE)


bestfit <-fit$best.model #Pick the best model using $
bestfit

fit

qqghyp(bestfit, gaussian = FALSE, line = TRUE, main=" ",xlab = " ",
       plot.legend = FALSE)

qqnorm(y)
plot(y)


set.seed(123)
x2 <- rt(250000, df = 9)
fitdistr(x2, "t", df = 9)
## allow df to vary: not a very good idea!
fitdistr(x2, "t")

v <- numeric(100)
for (i in 1:100)
{
  v[i] <- ks.test(x2, "pt", i)$p.value
}
plot(v, type="l")
which(max(v) == v)

