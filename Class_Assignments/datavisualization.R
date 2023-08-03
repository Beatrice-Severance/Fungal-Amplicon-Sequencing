data("mtcars")

plot(x = mtcars$wt, y = mtcars$mpg)

plot(x = mtcars$wt, y = mtcars$mpg, 
     ylab = "Miles per gallon", xlab = "Weight",
     font.lab = 6,
     pch = 20)

library(ggplot2)
#point adds points, smooth adds trend line; lm means linear model, se turns on/off confidence intervals
#switching the order can change how the data is visualized! if you move geom_point to the end, points will be on top
#scale_y_log10 changes the axis to log; depending on your data you can use this option
ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_smooth(method = lm, se = FALSE, color = "orange") +
  geom_point(aes(color = wt)) +
  xlab("Weight") +
  ylab("Miles per gallon") +
  theme_classic() +
  scale_color_gradient(low = "forestgreen", high = "limegreen")

#distributions can be best represented by box plots (according to Dr. Noel)
filename = "Bull_richness.csv"
bull.richness <- read.csv(filename, header = TRUE, na.strings = NA)

#subset your data to only include what you want to study
#The comma at the end is important because it states to select all columns! Without this there will be an error
bull.richness.soy.no.till <- bull.richness[bull.richness$Crop == "Soy" &
                                             bull.richness$Treatment == "No-till",]

#jitter puts dots on top of your box plot/bar chart in order to show a more full distribution
ggplot(bull.richness.soy.no.till, aes(x = GrowthStage, y = richness, color = Fungicide)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.9)) +
  xlab("Growth Stage") +
  ylab("Richness")

#use stat_summary to calculate statistics without having to manually calculate first
#position = "dodge" can make the bars not overlap on each other
#21 is non-filled circle
ggplot(bull.richness.soy.no.till, aes(x = GrowthStage, y = richness, fill = Fungicide)) +
  stat_summary(fun=mean, geom = "bar", position = "dodge") +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge") +
  xlab("Growth Stage") +
  ylab("Richness") +
  geom_point(position = position_jitterdodge(dodge.width = 0.9), shape = 21, color = "black") +
  scale_color_manual(values = c("blue", "green")) +
  scale_fill_manual(values = c("blue", "green"))

#\n is a line break; you can put this on an axis so that it splits up a long name if you'd like
ggplot(bull.richness.soy.no.till, aes(x = GrowthStage, y = richness, group = Fungicide, color = Fungicide)) +
  stat_summary(fun=mean, geom = "point") +
  stat_summary(fun=mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
  ylab("Richness")

# facet_wrap
#below code has error; will look into a fix
#bull.richness$GrowthStage <- factor(bull.richness$GrowthStage, levels = c("R3", "R4", "R6", "V6", "V8", "V15"))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(bull.richness, aes(x = GrowthStage, y = richness, group = Fungicide, color = Fungicide)) +
  stat_summary(fun=mean, geom = "point") +
  stat_summary(fun=mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
  ylab("Richness") +
  facet_wrap(~Crop*Treatment, scales = "free") +
  theme_classic() +
  scale_color_manual(values = cbbPalette) +
  geom_jitter(alpha = 0.5)

