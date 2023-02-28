#Setup
filename <- "MycotoxinData.csv"
mycotoxin <- read.csv(filename, na.strings = "na", header = TRUE)

library(ggplot2)
library(tidyverse)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Question 2
ggplot(mycotoxin, aes(x = Treatment, y = DON, fill = Cultivar)) +
  geom_boxplot() +
  xlab("") +
  ylab("DON (ppm)")

#Question 3
ggplot(mycotoxin, aes(x = Treatment, y = DON, fill = Cultivar)) +
  stat_summary(fun=mean, geom = "bar", position = "dodge") +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", width = 0.5) +
  xlab("") +
  ylab("DON (ppm)")

#Question 4
ggplot(mycotoxin, aes(x = Treatment, y = DON, fill = Cultivar)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.9), shape = 21, color = "black") +
  xlab("") +
  ylab("DON (ppm)")

ggplot(mycotoxin, aes(x = Treatment, y = DON, fill = Cultivar)) +
  stat_summary(fun=mean, geom = "bar", position = "dodge") +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", width = 0.5) +
  geom_point(position = position_jitterdodge(dodge.width = 0.9), shape = 21, color = "black") +
  xlab("") +
  ylab("DON (ppm)")

#Question 5
ggplot(mycotoxin, aes(x = Treatment, y = DON, fill = Cultivar)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.9), shape = 21, color = "black") +
  scale_fill_manual(values = c("#56B4E9", "#009E73")) +
  xlab("") +
  ylab("DON (ppm)")

ggplot(mycotoxin, aes(x = Treatment, y = DON, fill = Cultivar)) +
  stat_summary(fun=mean, geom = "bar", position = "dodge") +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", width = 0.5) +
  geom_point(position = position_jitterdodge(dodge.width = 0.9), shape = 21, color = "black") +
  scale_fill_manual(values = c("#56B4E9", "#009E73")) +
  xlab("") +
  ylab("DON (ppm)")

#Question 6
ggplot(mycotoxin, aes(x = Treatment, y = DON, fill = Cultivar)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.9), shape = 21, color = "black") +
  scale_fill_manual(values = c("#56B4E9", "#009E73")) +
  xlab("") +
  ylab("DON (ppm)") +
  facet_wrap(~Cultivar, scales = "free")

ggplot(mycotoxin, aes(x = Treatment, y = DON, fill = Cultivar)) +
  stat_summary(fun=mean, geom = "bar", position = "dodge") +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", width = 0.5) +
  geom_point(position = position_jitterdodge(dodge.width = 0.9), shape = 21, color = "black") +
  scale_fill_manual(values = c("#56B4E9", "#009E73")) +
  xlab("") +
  ylab("DON (ppm)") +
  facet_wrap(~Cultivar, scales = "free")

#Question 7
ggplot(mycotoxin, aes(x = Treatment, y = DON, fill = Cultivar)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.9), shape = 21, color = "black") +
  scale_fill_manual(values = c("#56B4E9", "#009E73")) +
  xlab("") +
  ylab("DON (ppm)") +
  facet_wrap(~Cultivar, scales = "free") +
  theme_minimal()

ggplot(mycotoxin, aes(x = Treatment, y = DON, fill = Cultivar)) +
  stat_summary(fun=mean, geom = "bar", position = "dodge") +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", width = 0.5) +
  geom_point(position = position_jitterdodge(dodge.width = 0.9), shape = 21, color = "black") +
  scale_fill_manual(values = c("#56B4E9", "#009E73")) +
  xlab("") +
  ylab("DON (ppm)") +
  facet_wrap(~Cultivar, scales = "free") +
  theme_classic()

#Question 8
ggplot(mycotoxin, aes(x = Treatment, y = DON, fill = Cultivar, alpha = 0.1)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.9), shape = 21, color = "black") +
  scale_fill_manual(values = c("#56B4E9", "#009E73")) +
  xlab("") +
  ylab("DON (ppm)") +
  facet_wrap(~Cultivar, scales = "free") +
  theme_minimal()

ggplot(mycotoxin, aes(x = Treatment, y = DON, fill = Cultivar, alpha = 0.1)) +
  stat_summary(fun=mean, geom = "bar", position = "dodge") +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", width = 0.5) +
  geom_point(position = position_jitterdodge(dodge.width = 0.9), shape = 21, color = "black") +
  scale_fill_manual(values = c("#56B4E9", "#009E73")) +
  xlab("") +
  ylab("DON (ppm)") +
  facet_wrap(~Cultivar, scales = "free") +
  theme_classic()

#Question 10
ggplot(mycotoxin, aes(x = Treatment, y = DON, fill = Cultivar, alpha = 0.1)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.9), shape = 21, color = "black") +
  scale_fill_manual(values = c("#56B4E9", "#009E73")) +
  xlab("") +
  ylab("DON (ppm)") +
  facet_wrap(~Cultivar*BioRep, scales = "free") +
  theme_classic()

ggplot(mycotoxin, aes(x = Treatment, y = DON, fill = Cultivar, alpha = 0.1)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.9), shape = 21, color = "black") +
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  xlab("") +
  ylab("DON (ppm)") +
  facet_wrap(~Treatment, scales = "free") +
  theme_classic()
