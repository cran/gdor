library(gdor)

data <- read.table(url("http://www.stat.umn.edu/geyer/gdor/catrec.txt"),
                  header = TRUE)

try(m1 <- gdor(y ~ (v1 + v2 + v3 + v4 + v5 + v6 + v7)^2, family=poisson,
               data=data))
try(m2 <- gdor(y ~ (v1 + v2 + v3 + v4 + v5 + v6 + v7)^3, family=poisson,
               data=data))
try(m3 <- gdor(y ~ (v1 + v2 + v3 + v4 + v5 + v6 + v7)^4, family=poisson,
               data=data))

try(anova(m1,m2))

try(anova(m2,m3))





