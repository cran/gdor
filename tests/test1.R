library(gdor)

x <- 1:30

y <- c(rep(0, 12), rep(1, 11), rep(0, 7)) 

m <- gdor(y ~ x + I(x^2), family=binomial)

summary(m)


