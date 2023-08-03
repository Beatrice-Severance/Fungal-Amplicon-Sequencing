z <- 1:200

mean(z)
sd(z)

vec <- z > 1
vec

df <- data.frame(z, vec)

colnames(df) <- c("z", "zlog")

df$zsquared <- z^2

dfsub <- subset(df, zsquared > 10 & zsquared < 100)
dfsub