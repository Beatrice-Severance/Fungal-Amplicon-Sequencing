#dependencies
library(tidyverse)
#install.packages("drc")
library(drc)

### Functions
#write a function to convert Fahrenheit to Celsius
(5*(degree_F - 32)/9)

#example
#function.name <- function(x) {
  # code goes here #
  # return(# output #) #
#}

FtoC <- function(x) {
  C <- (5*(x - 32)/9)
  return(C)
}

FtoC(80)

#don't use F as a part of the function because it is a logical operator (means FALSE!); use lowercase or something entirely different
CtoF <- function(x) {
  f <- (x*(9/5) + 32)
  return(f)
}

CtoF(0)

#for loops
for (i in 1:10) {
  print(i*2)
}

#example
for (i in -100:100) {
  result <- FtoC(i)
  print(result)
}

#try on a loaded dataset!
filename = "EC50.all.csv"
EC50.data <- read.csv(filename, header = TRUE)

#make a list of the data to loop
nm <- unique(EC50.data$is)

#seq_along() gives each observation a number to loop by
EC50.ll4 <- NULL
for (i in seq_along(nm)) {
  isolate1 <- drm(100 * EC50.data$relgrowth[EC50.data$is == nm[[i]]] ~
                    EC50.data$conc[EC50.data$is == nm[[i]]],
                  fct = LL.4(fixed = c(NA, NA, NA, NA),
                             names = c("Slope", "Lower", "Upper", "EC50")),
                  na.action = na.omit)
  summary.fit <- data.frame(summary(isolate1)[[3]])
  EC50 <- ED(isolate1, respLev = c(50), type = "relative",
             interval = "delta"[[1]])
  EC50
  isolate_ec_i <- data.frame(nm[i], EC50)
  colnames(isolate_ec_i) <- c("Isolate", "EC50")
  EC50.ll4 <- rbind.data.frame(EC50.ll4, isolate_ec_i)
}


ec50results <- EC50.data %>%
  group_by(is) %>%
  nest(ll.4.mod = map(data, ~drm(.$relgrowth ~ .$conc, 
                                 fct = LL.4()))) %>%
  mutate(ec50 = map(ll.4.mod, ~ED(., respLev = 50)[1])) %>%
  unnest(ec50) %>%
  ggplot(aes(ec50)) +
  geom_histogram()


ec50results$data[1]
ec50results$ll.4.mod[1]

  
