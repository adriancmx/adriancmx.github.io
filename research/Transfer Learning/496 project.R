library(KernSmooth)
library(locfit)
# Target domain Q
set.seed(123)
nq <- 200
x.q1 <- seq(0, 1, l = nq/2)
z.q1 <- rnorm(nq/2, 0, 1/3)
f <- function(x) sin(10*pi*x) + x^1.5 -0.1*x + max(0, 0.1- abs(x-0.5))
y.q1 <- f(x.q1) + z.q1
#y.q <- f(x.q1)
plot(x.q1, y.q1)
#lines(x.q1,y.q)


x.q2 <- seq(0, 1, l = nq/2)
z.q2 <- rnorm(nq/2, 0, 1/3)
y.q2 <- f(x.q2) + z.q2
plot(x.q2, y.q2)

# Source domain P
np <- 4800
x.p <- seq(0, 1, by = 1/np)
z.p <- rnorm(np, 0, 1/3 )
g <-  function(x) sin(10*pi*x) + x^1.5 -0.1*x
y.p <- g(x.p) +  z.p
plot(x.p, y.p)

# estimate f.ref,
f.hat.ref <- locfit(y.q1~lp(x.q1, deg=3, h=1/trunc((nq/2)^(1/(2*1+1)))),kern='gauss')
plot(f.hat.ref,band = 'global',col = 'red ')
points(x.q1, y.q1, col = 'black')
pred.ref <- predict(f.hat.ref, x.q1, se = TRUE, interval = "confidence", level = 0.95)

pred.ref2 <- predict(f.hat.ref, x.q2, se = TRUE, interval = "confidence", level = 0.95)
LCI <- 2*qnorm(0.975)*pred.ref$se.fit
##lower <- pred$fit - qnorm(0.975)*pred$se.fit
##upper <- pred$fit + qnorm(0.975)*pred$se.fit
##lines(x.q1,lower, col = 'green')
##lines(x.q1,upper, col = 'blue')

mse.lpr <- 1/(nq/2)*sum((y.q1-pred.ref$fit)^2)
mse.lpr

# estimate f.raw
f.hat.raw <- locfit(y.p~lp(x.p, deg=3, h=1/trunc(np^(1/(2*1+1)))),kern='gauss')
pred.raw  <- predict(f.hat.raw, x.q1, se = TRUE, interval = "confidence", level = 0.95)
pred.raw2 <- predict(f.hat.raw, x.q2, se = TRUE, interval = "confidence", level = 0.95)
plot(f.hat.raw, band = 'global',col = 'red ')
points(x.p, y.p, col = 'black')

# Define the objective function
obj <- function(psi, B, x, y,pred.raw, pred.ref, LCI) {
  y_hat <- pred.ref + sign(pred.raw + B %*% psi- pred.ref ) * pmin(abs(pred.raw + B %*% psi- pred.ref), LCI/2)
  sum((y - y_hat)^2)
}

# Set initial value for psi
B <- poly(x.q2, 3)
psi_init <- c(1,1,1)

result <- optim(psi_init, obj, B = B, x=x.q2,y=y.q2, pred.raw= pred.raw2$fit, pred.ref=pred.ref2$fit, LCI=LCI)

# Extract the best psi
psi_best <- result$par

# Evaluate the fitted values with the best psi
f.hat.ct <- pred.ref$fit + sign(pred.raw$fit + poly(x.q1,3)%*%psi_best - pred.ref$fit ) * 
  pmin(abs(pred.raw$fit + poly(x.q1,3)%*%psi_best - pred.ref$fit), LCI/2)

mse.act <- 1/(nq/2)*sum((y.q1-f.hat.ct)^2)
mse.act

plot(x.q1,y.q1)
lines(x.q1,f.hat.ct)
