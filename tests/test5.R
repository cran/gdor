library(gdor) 
team.names <- c("ants", "beetles", "cows", "dogs",
"egrets", "foxes", "gerbils", "hogs")

data <- matrix(c(NA, 2, 2, 2, 2, 2, 2, 2, 0, NA,
1, 2, 2, 2, 2, 2, 0, 1, NA, 2, 1, 2, 2, 2, 0,
0, 0, NA, 1, 1, 2, 2, 0, 0, 1, 1, NA, 1, 2, 2,
0, 0, 0, 1, 1, NA, 2, 2, 0, 0, 0, 0, 0, 0, NA,
1, 0, 0, 0, 0, 0, 0, 1, NA), byrow = TRUE, nrow = 8)
dimnames(data) <- list(team.names, team.names)

wins <- data[upper.tri(data)]
team.plus <- row(data)[upper.tri(data)]
team.minus <- col(data)[upper.tri(data)]
modmat <- matrix(0, length(wins), nrow(data))

for (i in 1:ncol(modmat)) {
modmat[team.plus == i, i] <- 1
modmat[team.minus == i, i] <- (-1)
}
losses <- 2 - wins
resp <- cbind(wins, losses)

# out <- glm(resp ~ modmat + 0, family = binomial)
out <- gdor(resp ~ modmat + 0, family = binomial) 

summary(out)