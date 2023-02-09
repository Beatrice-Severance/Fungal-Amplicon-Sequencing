#R is a fungi

2+2
2-2
2/2
2*2

x = 3
y <- 2

x+y
2 -> z

name <- "Beatrice"
seven <- "7"

class(seven)
class(x)

vec <- c(1,2,3,4,5,6,7) #numeric vector
vec = c(1:7)
vec <- 1:7

vec2 <- c("Zach", "Jie", "Mitch")
vec3 <- c(TRUE, FALSE, TRUE)

vec2[2]
vec + x

mean(vec)
sd(vec)
sd(vec)/sqrt(7) #standard error

t <- 1:10

t[(t > 8) | (t < 5)]
t[t != 2]

1 %in% t

mat1 <- matrix(data = c(1, 2, 3), nrow = 3, ncol = 3)

class(mat1)

mat2 <- matrix(data = c("Zach", "Jie", "Mitch"), nrow = 3, ncol = 3)

mat1[5]
mat1[2, 2] #rows comma columns
mat2[,1]

df <- data.frame(mat1[,1], mat2[,1])
colnames(df) <- c("value", "name")

df[,1]
df$value

df[df$name == "Jie",] #all rows associated with the specified name
df$value[df$name == "Jie"] #value associated with the name Jie

subset(df, name == "Jie")
