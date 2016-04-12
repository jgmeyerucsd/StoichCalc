require(graphics)

dnorm(0) == 1/sqrt(2*pi)
dnorm(1) == exp(-1/2)/sqrt(2*pi)
dnorm(1) == 1/sqrt(2*pi*exp(1))


par(mfrow = c(2,1))
plot(function(x) dnorm(x, log = TRUE), -60, 50,
     main = "log { Normal density }")
curve(log(dnorm(x)), add = TRUE, col = "red", lwd = 2)

x   <- seq(5,15,length=1000)
y   <- dnorm(x,mean=10, sd=1.3)
y   <- 0.6*dnorm(x,mean=10, sd=1)
par(lwd=10)
plot(x,y, type="l", lwd=5, col="black",xlim=c(0,20))
plot(x,y, type="l", lwd=1,xlim=c(0,20))

lines(x,0.8*y,xlim=c(0,20),lwd=5)
lines(x,0.6*y, type="l", lwd=1,xlim=c(0,20))


dev.off()
