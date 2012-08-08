library(gdor)

m1 <- gdor(cbind(ncases, ncontrols) ~ agegp + tobgp * alcgp, 
	   family = binomial, data=esoph)

summary(m1)
