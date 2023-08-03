#load libraries
library(tidyverse)
library(lme4)
library(emmeans)
#install.packages("multcomp")
library(multcomp)
#install.packages("multcompView")
library(multcompView)

#practice
data("mtcars")

ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_smooth(method = lm, se = FALSE) +
  geom_point()

#what are the stats behind lm?
#low p-value means that there is a significant relationship between the two variables
lm1 <- lm(mpg~wt, data = mtcars)
summary(lm1)
#make an anova table to show the results; the p-values are similar! an ANOVA is a regression and vice versa
anova(lm1)
#R-squared and adjusted R-squared: percent variation in y that is explained by the variation in x

#run a correlation test; wow the p-value is still the same!
cor.test(mtcars$wt, mtcars$mpg)

#evaluate residuals; check assumptions; maybe the relationship isn't completely linear but it's still fair
ggplot(lm1, aes(y = .resid, x = .fitted)) +
  geom_point() +
  geom_hline(yintercept = 0)

#the below command shows multiple plots so you can get a quick view of your data
plot(lm1)

#load bull richness
filename = "./Class Assignments/Bull_richness.csv"
bull.rich <- read.csv(filename, header = TRUE)

bull.rich %>%
  filter(GrowthStage == "V8" & Treatment == "Conv.") %>%
  ggplot(aes(x = Fungicide, y = richness)) +
  geom_boxplot()

bull.rich.sub <- bull.rich %>%
  filter(GrowthStage == "V8" & Treatment == "Conv.")

t.test(richness~Fungicide, data = bull.rich.sub, var.equal = TRUE)

summary(lm(richness~Fungicide, data = bull.rich.sub))

anova(lm(richness~Fungicide, data = bull.rich.sub))

bull.rich.sub2 <- bull.rich %>%
  filter(Fungicide == "C" & Treatment == "Conv." & Crop == "Corn")

ggplot(bull.rich.sub2, aes(x = GrowthStage, y = richness)) +
  geom_boxplot()

lm.growth <- lm(richness~GrowthStage, data = bull.rich.sub2)
summary(lm.growth)

anova(lm.growth)

#choose what you would like to separate by with the ~
#multcomp:: add this if libraries are conflicting for some reason
#emmeans is estimated marginal means
lsmeans <- emmeans(lm.growth, ~GrowthStage)
results <- cld(lsmeans, alpha = 0.05, reversed = TRUE, details = TRUE)
results

### Interactions

bull.rich.sub3 <- bull.rich %>%
  filter(Treatment == "Conv." & Crop == "Corn")

#long way to write
lm.inter <- lm(richness ~ GrowthStage + Fungicide + GrowthStage:Fungicide, data = bull.rich.sub3)
#short way! test each individual term plus the interaction term
lm.inter <- lm(richness ~ GrowthStage*Fungicide, data = bull.rich.sub3)
summary(lm.inter)
anova(lm.inter)

ggplot(bull.rich.sub3, aes(x = GrowthStage, y = richness, fill = Fungicide)) +
  geom_boxplot()

#| means within
lsmeans <- emmeans(lm.inter, ~Fungicide|GrowthStage)
results <- cld(lsmeans, alpha = 0.05, reversed = TRUE, details = TRUE)
results

#what if you want it the other way around?
lsmeans <- emmeans(lm.inter, ~GrowthStage|Fungicide)
results <- cld(lsmeans, alpha = 0.05, reversed = TRUE, details = TRUE)
results


lm.inter <- lm(richness ~ GrowthStage*Fungicide, data = bull.rich.sub3)
#Add a random effect with the +, the lmer function needs this (show mixed effects)
lme1 <- lmer(richness ~ GrowthStage*Fungicide + (1|Rep), data = bull.rich.sub3)
summary(lm.inter)
summary(lme1)
#standard errors are going down a bit, so the estimates are more confident because we included the random factor

