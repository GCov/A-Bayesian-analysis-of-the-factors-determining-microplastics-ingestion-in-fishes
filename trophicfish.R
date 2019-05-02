trophicfish <- 
  read.csv("~/GRAD school/Trophic Fish Paper/Data/trophicfishupdate.csv")

summary(trophicfish)

trophicfish$year <- as.factor(trophicfish$year)

# concentration by weight

trophicfish$conc <- trophicfish$Mpsgut/trophicfish$W
summary(trophicfish)
trophicfish$IR <- trophicfish$ingest.rate/100

library(plyr)

## Consolidate averages by study

trophicfish2 <- ddply(trophicfish, 
                      c('author', 'year', 'region', 'species', 'family',
                        'common.name', 'genus', 'environment', 'climate', 
                        'red.list', 'feeding.type', 'feeding.habit', 'TL', 
                        'min.size', 'float.meth', 'dig.meth', 'count.meth', 
                        'polymer.meth', 'maj.fib', 'maj.under.one.mm', 
                        'maj.col', 'exclude.fib'), 
                      summarise, 
                      fork.length = sum(fork.length*N)/sum(N), 
                      total.length = sum(total.length*N)/sum(N), 
                      W = sum(W*N)/sum(N), 
                      GW = sum(GW*N)/sum(N),
                      N = sum(N), 
                      IR = sum(IR*N)/sum(N), 
                      Mpsgut = sum(Mpsgut*N)/sum(N), 
                      SDMPsgut = sum(SDMPsgut*N)/sum(N)
                      )

trophicfish2$IR <- with(trophicfish2, ingested/N)
summary(trophicfish2)
length(trophicfish2$species) # 482 dta points
length(trophicfish$species) # consolidated from 613 data point
length(unique(trophicfish2$species)) # 383 species
length(unique(trophicfish2$family)) # from 131 families

trophicfish2$study <- with(trophicfish2, paste0(author, year))
trophicfish2$study <- as.factor(trophicfish2$study)
head(trophicfish2)
length(unique(trophicfish2$study))  # 63 studies

summary(trophicfish2$min.size)

trophicfish2$min.size <- mapvalues(trophicfish2$min.size, from = c(""), 
                                   to = c(NA))

summary(trophicfish2$min.size)

trophicfish2$min.size <- 
  as.numeric(levels(trophicfish2$min.size))[trophicfish2$min.size]

summary(trophicfish2$min.size)
hist(trophicfish2$min.size)

trophicfish2$climate <- 
  mapvalues(trophicfish2$climate, from = "Benthopelagic", to = "Temperate")

gutdata <- subset(trophicfish2, Mpsgut != 'NA')
summary(gutdata)
summary(gutdata$author)
length(gutdata$species) # 366 data points
length(unique(gutdata$species)) # 382 species
length(unique(gutdata$family)) # from 120 families
length(unique(gutdata$study)) # from 49 studies

summary(gutdata)
gutdata$region <- as.character(gutdata$region)
gutdata$region <- as.factor(gutdata$region)
summary(gutdata$region)

gutdata$environment <- as.character(gutdata$environment)
gutdata$environment <- as.factor(gutdata$environment)
summary(gutdata$environment)

summary(gutdata$TL)
gutdata$TL <- as.numeric(gutdata$TL)

summary(gutdata$Mpsgut)
gutdata$Mpsgut <- as.numeric(gutdata$Mpsgut)

library(ggplot2)

ggplot(gutdata, aes(x=climate , y=log(Mpsgut))) + 
  geom_boxplot(size=1, fill = 'grey80') + 
  facet_wrap(gutdata$environment) +
  geom_smooth(method = 'lm', colour = "black", se = TRUE) +
  xlab('Trophic Level') + 
  ylab("# of MPs in gut") +
  theme_bw() +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(size = 16, angle = 50, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

gutdata$gutconc <- with(gutdata, Mpsgut/log(W))
summary(gutdata$gutconc)

gut.conc <- subset(gutdata, gutconc != 'NA' & gutconc != 'Inf')
summary(gut.conc)

length(gut.conc$species) # 191 data points
length(unique(gut.conc$species)) # 169 species
length(unique(gut.conc$family)) # from 70 families
length(unique(gut.conc$study)) # 24 studies
length(unique(gut.conc$region)) ## 13 regions

summary(gut.conc)

gutdata$feeding.habit <- mapvalues(gutdata$feeding.habit, 
                                   from = levels(gutdata$feeding.habit),
                                   to = c("Not listed", 
                                          "Browsing on substrate", 
                                          "Filtering plankton",
                                          "Grazing on aquatic plants", 
                                          "Hunting macrofauna",
                                          "Selective plankton feeding", 
                                          "Variable"))

levels(gutdata$feeding.habit)

gutdata$feeding.habit <- factor(
  gutdata$feeding.habit,
  levels = c(
    "Not listed",
    "Browsing on substrate",
    "Grazing on aquatic plants",
    "Filtering plankton",
    "Selective plankton feeding",
    "Hunting macrofauna",
    "Variable"
  )
)
levels(gutdata$feeding.habit)

hist(gutdata$Mpsgut)

plot(log(Mpsgut + 1) ~ N, data = gutdata)  # variance decreases with N

gutdata$negN <- 0 - gutdata$N

plot(log(Mpsgut + 1) ~ negN, data = gutdata)  # seems to solve the issue

plot(log(Mpsgut + 1) ~ TL, data = gutdata) # variance doesn't look too bad

library(nlme)
library(MuMIn)

## First make trophic level/allometric models

M1 <- lme(log(Mpsgut+1) ~ TL, random = ~1|study, 
          weights = varPower(form =~ N), 
          data = gutdata, method = "REML")
summary(M1)
plot(resid(M1) ~ fitted(M1))
abline(0,0)

## Compare fit without weights

M2 <- lme(log(Mpsgut+1) ~ TL, random = ~1|author/year, 
          data = gutdata, method = "REML")

AICc(M1, M2)  # better fit with weights

plot(resid(M2) ~ fitted(M2))

plot(resid(M1) ~ gutdata$TL)
plot(resid(M1) ~ log(gutdata$N))

plot(Mpsgut ~ TL, data = gutdata)
lines(exp(predict(M1, type = 'response')) - 1, col = 'red')

summary(M1)  # no correlation between trophic level and MPS, p = 0.27

# switch to ML

M3 <- lme(log(Mpsgut+1) ~ TL, random = ~1|study, 
          weights = varPower(form =~ N), 
          data = gutdata, method = "ML")

summary(M3)  # same p-value

## Plot

png(
  filename = "trophiclevelplot.png",
  width = 9,
  height = 8,
  units = "cm",
  pointsize = 7,
  res = 600
)

ggplot(gutdata, aes(x=TL , y=Mpsgut)) +
  geom_jitter(aes(size=gutdata$N), alpha = 5/10, shape = 21) + 
  xlab("") +
  scale_x_continuous(limits = c(1.9, 4.75), 
                     breaks = c(2, 2.5, 3, 3.5, 4, 4.5, 5)) +
  ylab("# MPs in Gut") +
  xlab("Trophic Level") +
  scale_size_continuous(range = c(0.1, 4)) +
  scale_y_continuous(trans = 'log1p',
                     breaks = c(0, 1, 10, 20, 30)) +
  theme_classic() +
  labs(size="n") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme(text = element_text(size=7), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="right",
        legend.justification="left",
        legend.margin=margin(0,0,0,0))

dev.off()

## Now consider effect of allometry

allo <- subset(gutdata, total.length != 'NA')

length(allo$total.length)  # 172 data point remaining

M4 <- lme(log(Mpsgut+1) ~ TL*total.length, random = ~1|study, 
          weights = varPower(form =~ N), 
          data = allo, method = "REML")

summary(M4)
plot(resid(M4) ~ fitted(M4))  # variance does not look homogenous

## try without variance structure

M5 <- lme(log(Mpsgut+1) ~ TL*total.length, random = ~1|study,
          data = allo, method = "REML")

plot(resid(M5) ~ fitted(M5))  # doesn't look much better

AICc(M4, M5)  # fit is a little better with error structure

# figure out what is driving the increasing variance

plot(resid(M4) ~ allo$TL)
plot(resid(M4) ~ allo$total.length)  # looks like it's total length

# try log-transformaing total length

M6 <- lme(log(Mpsgut+1) ~ TL*log(total.length), random = ~1|study, 
          weights = varPower(form =~ N), 
          data = allo, method = "REML")

plot(resid(M6) ~ fitted(M6))  # looks better

AICc(M4, M6)  # M6 is a much better fit

summary(M6)  # nothing is significant

plot(exp(predict(M6)) ~ allo$TL)
plot(exp(predict(M6)) ~ allo$total.length)

## refit with ML

M7 <- lme(log(Mpsgut+1) ~ TL*log(total.length), random = ~1|study, 
          weights = varPower(form =~ N), 
          data = allo, method = "ML")

summary(M7)

png(
  filename = "totallengthplot.png",
  width = 9,
  height = 8,
  units = "cm",
  pointsize = 7,
  res = 600
)

ggplot(allo, aes(x=total.length , y=Mpsgut)) +
  geom_jitter(aes(size=N), alpha = 5/10, shape = 21) + 
  xlab("") +
  ylab("# MPs in Gut") +
  xlab("Total Length (cm)") +
  scale_size_continuous(range = c(1, 4)) +
  scale_x_continuous(trans = 'log10',
                     breaks = c(0, 1, 10, 100, 1000),
                     limits = c(0.1, 1000)) +
  scale_y_continuous(trans = 'log1p',
                     breaks = c(0, 1, 5, 10, 15)) +
  theme_classic() +
  labs(size="n") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme(text = element_text(size=7), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="right",
        legend.justification="left",
        legend.margin=margin(0,0,0,0))

dev.off()

## Try all of this again in terms of concentration by weight

M8 <- lme(log(gutconc+1) ~ TL, random = ~1|study, 
          weights = varPower(form =~ N), 
          data = gut.conc, method = "REML")

plot(resid(M8) ~ fitted(M8))
plot(resid(M8) ~ gut.conc$TL)
plot(resid(M8) ~ gut.conc$N)

## Try without weights

M9 <- lme(log(gutconc+1) ~ TL, random = ~1|study, 
          data = gut.conc, method = "REML")
AICc(M8, M9)  # fits better with weights

M10 <- lme(log(gutconc+1) ~ TL, random = ~1|study, 
          weights = varPower(form =~ N), 
          data = gut.conc, method = "ML")

summary(M10)  # p = 0.07

## account for body size 

allo2 <- subset(gut.conc, total.length != 'NA')
length(allo2$total.length)  # 146 data points

M11 <- lme(log(gutconc+1) ~ TL*log(total.length), random = ~1|study, 
           weights = varPower(form =~ N), 
           data = allo2, method = "REML")

summary(M11)
plot(resid(M11) ~ fitted(M11))
plot(resid(M11) ~ allo2$TL)
plot(resid(M11) ~ log(allo2$total.length))
plot(resid(M11) ~ log(allo2$N))

M12 <- lme(log(gutconc+1) ~ TL*log(total.length), random = ~1|study, 
           weights = varPower(form =~ N), 
           data = allo2, method = "ML")

summary(M12)

png(
  filename = "concplot.png",
  width = 9,
  height = 8,
  units = "cm",
  pointsize = 7,
  res = 600
)

ggplot(allo2, aes(x=TL , y=gutconc)) +
  geom_jitter(aes(size=N), alpha = 5/10, shape = 21) + 
  xlab("") +
  ylab("# MPs in Gut per Gram Body Weight") +
  xlab("Trophic Level") +
  scale_y_continuous(trans = 'log1p',
                     breaks = c(0:6)) +
  theme_classic() +
  labs(size="n") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme(text = element_text(size=7), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="right",
        legend.justification="left",
        legend.margin=margin(0,0,0,0))

dev.off()

## Now look at habitat use and foraging method

summary(gutdata$environment)
gutdata$environment[gutdata$environment == 'Pelagic'] <- 'Pelagic-neritic'
gutdata$environment <- as.character(gutdata$environment)
gutdata$environment <- as.factor(gutdata$environment)

summary(gutdata$climate)

summary(gutdata$feeding.habit)

habit <- subset(gutdata, feeding.habit != 'Not listed')
summary(habit$feeding.habit)
habit$feeding.habit <- as.character(habit$feeding.habit)
habit$feeding.habit <- as.factor(habit$feeding.habit)

M13 <- lme(log(Mpsgut+1) ~ environment + climate + feeding.habit, 
           random = ~1|study, 
           weights = varPower(form =~ N), 
           data = habit, method = "REML")
plot(resid(M13) ~ fitted(M13))

summary(M13)

## refit with ML

M14 <- lme(log(Mpsgut+1) ~ environment + climate + feeding.habit, 
           random = ~1|study, 
           weights = varPower(form =~ N), 
           data = habit, method = "ML")

## test significance of the predicotrs

M15 <- lme(log(Mpsgut+1) ~ climate + feeding.habit, 
           random = ~1|study, 
           weights = varPower(form =~ N), 
           data = habit, method = "ML")

M16 <- lme(log(Mpsgut+1) ~ environment + feeding.habit, 
           random = ~1|study, 
           weights = varPower(form =~ N), 
           data = habit, method = "ML")

M17 <- lme(log(Mpsgut+1) ~ environment + climate, 
           random = ~1|study, 
           weights = varPower(form =~ N), 
           data = habit, method = "ML")

anova(M14, M15)  # environment not sig, p=0.76
anova(M14, M16)  # climate not sig. p=0.53
anova(M14, M17)  # feeding habit not sig. p = 0.64

## Now test effect of region

summary(trophicfish2$IR)

ingestion <- subset(trophicfish2, IR != "NA")

summary(ingestion$IR)

length(ingestion$species) # 413 data points
length(unique(ingestion$species)) # 334 species
length(unique(ingestion$family)) # from 128 families
length(unique(ingestion$study)) # 56 studies

levels(ingestion$region)

summary(ingestion$region)

ingestion$region[ingestion$region == 'Northeast Atlantic'] <-
  'Atlantic, Northeast'
ingestion$region <- as.character(ingestion$region)
ingestion$region <- as.factor(ingestion$region)

library(lme4)

summary(ingestion)

M18 <- glmer(IR ~ region + environment + TL +
               (1|study), 
           data = ingestion, 
           family = binomial)

## try removing variables

drop1(M18)  # remove environment

M19 <- glm(IR ~ region + TL, 
           data = ingestion, 
           family = binomial)

AICc(M18, M19)

drop1(M19)  # remove TL

M20 <- glm(IR ~ region, 
           data = ingestion, 
           family = binomial)

plot(resid(M20) ~ fitted(M20))
abline(0,0)

plot(resid(M20) ~ ingestion$region)

summary(M20)

ingestion$region <- relevel(ingestion$region, 'Pacific, Southwest')

png(
  filename = "Figure 1.png",
  width = 14,
  height = 13,
  units = "cm",
  pointsize = 7,
  res = 600
)

ggplot(ingestion) +
  geom_count(aes(x = reorder(region, IR, mean), 
                 y = predict(M20, type = 'response'))) +
  geom_errorbar(aes(
    x = reorder(region, IR, sum),
    ymin = 
      predict(M20, type = 'response') - 
      predict(M20, type = 'response',  se.fit = TRUE)$se.fit,
    ymax = 
      predict(M20, type = 'response') + 
      predict(M20, type = 'response', se.fit = TRUE)$se.fit),
    size = 0.5
    ) +
  labs(y = 'Ingestion Rate', x = 'FAO Region', size = 'Number of Datapoints') +
  coord_flip(ylim = c(-0.06,1.2)) +
  scale_y_continuous(breaks = c(0,0.25, 0.5, 0.75,1)) +
  theme_classic() + 
  theme(text = element_text(size=7), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = c(1, 1)) +
  annotate("text", x=14.4, y=0.871, label = "A", size = 3) +
  annotate("text", x=13.4, y=0.775, label = "A,B,C,D,E,F", size = 3) +
  annotate("text", x=12.4, y=0.59, label = "A,C,D,F", size = 3) +
  annotate("text", x=11.4, y=0.492, label = "C", size = 3) +
  annotate("text", x=10.4, y=0.41, label = "B,C,D", size = 3) +
  annotate("text", x=9.4, y=0.4, label = "C,D", size = 3) +
  annotate("text", x=8.5, y=0.294, label = "D", size = 3) +
  annotate("text", x=7.4, y=0.288, label = "B,C,D,E", size = 3) +
  annotate("text", x=6.4, y=0.214, label = "B,C,D,E", size = 3) +
  annotate("text", x=5.5, y=0.155, label = "B,E", size = 3) +
  annotate("text", x=4.5, y=0.154, label = "B,F", size = 3) +
  annotate("text", x=3.4, y=0.095, label = "B,C,D,E", size = 3) +
  annotate("text", x=2.5, y=0.068, label = "E,F", size = 3) +
  annotate("text", x=1.4, y=0.055, label = "E,F", size = 3)

dev.off()
#####






A <-
ggplot(gutdata, aes(x=TL , y=log(Mpsgut+1))) + 
  geom_smooth(aes(x=gutdata$TL, y = M1predict), 
            col = "red", size = 1, alpha = 5/10) +
  geom_point(aes(size=gutdata$N), alpha = 5/10) + 
  xlab("") +
  scale_x_continuous(limits = c(1.9, 4.75), breaks = c(2, 2.5, 3, 3.5, 4, 4.5, 5)) +
  ylab("ln(MPs + 1)") +
  theme_bw() +
  labs(size="n") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme(text = element_text(size=12), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="right",
        legend.justification="left",
        legend.margin=margin(0,0,0,0))


B <-
ggplot(gutdata, aes(x = TL , y = log(Mpsgut+1), colour = feeding.habit)) +
  geom_jitter(aes(size=gutdata$N), alpha = 5/10) + 
  geom_line(aes(x = TL, y = predict.av, colour = feeding.habit, alpha = 5/10), size = 0.8) +
  xlab('Trophic Level') +
  scale_x_continuous(limits = c(1.9, 4.75), breaks = c(2, 2.5, 3, 3.5, 4, 4.5, 5)) +
  ylab("ln(MPs + 1)") +
  theme_bw() +
  scale_alpha(guide = FALSE) +
  scale_size_area(guide = FALSE) +
  scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9", 
                              "#009E73", "#CC79A7", "#0072B2", "#D55E00")) +
  labs(size="n", colour="Feeding Habit") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme(text = element_text(size=12), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="right",
        legend.justification="left",
        legend.margin=margin(0,0,0,0))

library(cowplot)

png(filename = "gutMPsplots.png",width = 24, height = 20, units = "cm", pointsize = 12, res=600)

plot_grid(A, B, labels=c("A","B"), ncol = 1, hjust = -1.4, nrow =2, label_size=12, align = 'hv')

dev.off()



# plot by weight

levels(gut.conc$feeding.habit)

gut.conc$feeding.habit <- mapvalues(gut.conc$feeding.habit, from = levels(gut.conc$feeding.habit),
                                   to = c("Not listed", "Browsing on substrate", "Filtering plankton",
                                          "Grazing on aquatic plants", "Hunting macrofauna",
                                          "Selective plankton feeding", "Variable"))

levels(gut.conc$feeding.habit)

gut.conc$feeding.habit <- factor(gut.conc$feeding.habit, levels = c("Not listed", "Browsing on substrate", 
                                                                  "Grazing on aquatic plants",
                                                                  "Filtering plankton", 
                                                                  "Selective plankton feeding", 
                                                                  "Hunting macrofauna",
                                                                  "Variable"))

png(filename = "weightplot.png",width = 24, height = 20, units = "cm", pointsize = 12, res=600)

ggplot(gut.conc, aes(x=TL , y=log(W), colour = feeding.habit)) + 
  geom_line(aes(x=TL, y = fitted(BSM.2)), col = "black", size = 1, alpha = 3/10) + 
  geom_line(aes(x=TL, y = (0.8762 + (1.5777*TL))), size = 1, linetype = 'dashed', 
            col = 'grey50', alpha = 3/10) +
  geom_line(aes(x=TL, y = (-1.4262 + (0.9211*TL))), size = 1, linetype = 'dashed', 
            col = 'grey50', alpha = 3/10) +
  geom_jitter(aes(size=gut.conc$N), alpha = 5/10) + 
  xlab('Trophic Level') + 
  ylab("ln(Weight) (g)") +
  theme_bw() +
  scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9", 
                              "#009E73", "#CC79A7", "#0072B2", "#D55E00")) +
  labs(size="n", colour="Feeding Habit") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  scale_x_continuous(limits = c(1.9, 4.75), breaks = c(2, 2.5, 3, 3.5, 4, 4.5, 5)) +
  theme(text = element_text(size=12), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

dev.off()


#####################################
summary(trophicfish2$IR)

ingestion <- subset(trophicfish2, IR != "NA")

summary(ingestion$IR)

length(ingestion$species) # 309 data points
length(unique(ingestion$species)) # 263 species
length(unique(ingestion$family)) # from 111 families
length(unique(ingestion$author)) # 39 studies
tapply(ingestion$year, ingestion$author, unique)
levels(ingestion$region)

summary(ingestion$region)

library(car)
library(lme4)
library(nlme)
library(mgcv)

summary(ingestion)

M2 <- gam(IR ~ region + environment + TL, data = ingestion, weights = N, family = "binomial")
summary(M2)
anova(M2)
plot(resid(M2) ~ fitted(M2))
plot(resid(M2) ~ ingestion$region)
plot(resid(M2) ~ ingestion$environment)
plot(resid(M2) ~ ingestion$TL)

plot(IR ~ TL, data = ingestion, ylim = c(0,1))
points(predict(M2, type = "response") ~ ingestion$TL, col = "red")
lines(spline(predict(M2, type = "response") ~ ingestion$TL), col = "red")


ingestion$region <- relevel(ingestion$region, "Western Central Pacific")


par(mfrow=c(1,2))
plot(IR ~ reorder(region, IR, mean), data = ingestion)
plot(predict(M2, method = "response") ~ reorder(ingestion$region, ingestion$IR, mean), col = "red")

levels(ingestion$region)
ingestion$region <- relevel(ingestion$region, "Gulf of Mexico")

tapply(ingestion$IR, ingestion$region, mean)
tapply(ingestion$IR, ingestion$region, sd)
tapply(ingestion$IR, ingestion$region, length)
tapply(ingestion$author, ingestion$region, unique)


png(filename = "ingestionrateplot.png",width = 25, height = 20, units = "cm", pointsize = 12, res=600)

ggplot(ingestion, aes(x=reorder(region, IR, mean) , y=IR)) + 
  geom_boxplot(size=1, varwidth = TRUE, fill = 'grey80') + 
  xlab('Region') + 
  ylab("MP Ingestion Rate") +
  theme_bw() +
  theme(text = element_text(size=12), 
        axis.text.x = element_text(size = 12, angle = 50, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  annotate("text", x=1, y=0.1, label = "A,B,C,D,E,F") +
  annotate("text", x=2, y=0.2, label = "A") +
  annotate("text", x=3, y=0.25, label = "A") +
  annotate("text", x=4, y=0.7, label = "A,B") +
  annotate("text", x=5, y=1.05, label = "B") +
  annotate("text", x=6, y=0.6, label = "B,C") +
  annotate("text", x=7, y=1.05, label = "A") +
  annotate("text", x=8, y=0.7, label = "C") +
  annotate("text", x=9, y=0.55, label = "C") +
  annotate("text", x=10, y=1.05, label = "D") +
  annotate("text", x=11, y=1.05, label = "E") +
  annotate("text", x=12, y=0.9, label = "E") +
  annotate("text", x=13, y=1.05, label = "F")
  
dev.off()


summary(ingestion$min.size)
summary(ingestion)
ggplot(ingestion, aes(x=min.size , y=IR)) + 
  geom_point(size=1) + 
  facet_wrap(ingestion$region) +
  geom_smooth(method = 'lm', colour = "black", se = TRUE) +
  xlab('Minimum particle size') + 
  ylab("Ingestion rate") +
  theme_bw() +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(size = 16, angle = 50, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

minsize <- gutdata
summary(minsize$min.size)
minsize$min.size[minsize$min.size > 0 & minsize$min.size < 100] <- "<100"
minsize$min.size[minsize$min.size >= 100 & minsize$min.size < 500] <- "100-500"
minsize$min.size[minsize$min.size >= 500 & minsize$min.size < 1000] <- "500-1000"
minsize$min.size[minsize$min.size >= 1000] <- ">1000"

minsize$min.size <- as.factor(minsize$min.size)
summary(minsize$min.size)
levels(minsize$min.size)
minsize$min.size <- mapvalues(minsize$min.size, from = "0", to = "No limit")
levels(minsize$min.size)

plot(Mpsgut ~ min.size, data = minsize)
plot(TL ~ min.size, data = minsize)

plot(min.size ~ TL, data = minsize)
plot(min.size ~ region, data = minsize)

SM.1 <- lmer(log(Mpsgut+1) ~ min.size + (1|author/year), data = minsize, weights = N)
plot(SM.1)
summary(SM.1)
anova(SM.1)
SM.2 <- lmer(log(Mpsgut+1) ~ (1|author/year), data = minsize, weights = N)
anova(SM.1, SM.2)

plot(maj.under.one.mm ~ TL, data = trophicfish2)
plot(maj.fib ~ TL, data = trophicfish2)
plot(Mpsgut ~ polymer.meth, data = trophicfish2)
plot(IR ~ polymer.meth, data = ingestion)
plot(Mpsgut ~ dig.meth, data = trophicfish2)
tapply(trophicfish2$Mpsgut, trophicfish2$dig.meth, length)
plot(dig.meth ~ region, data = trophicfish2)

write.csv(trophicfish2, "trophicfishprocessed.csv")

tapply(minsize$author, minsize$min.size, unique)

minsize$min.size <- factor(minsize$min.size, levels = c("No limit", "<100", "100-500", ">1000"))

png(filename = "minimumsizeeffect.png",width = 24, height = 20, units = "cm", pointsize = 12, res=600)

ggplot(minsize, aes(x=min.size , y=log(Mpsgut + 1))) + 
  geom_boxplot(size=1, fill = 'grey80') + 
  xlab(expression("Minimum"~"particle"~"detection"~"limit"~"("*mu*"m)")) + 
  ylab("ln(MPs + 1)") +
  theme_bw() +
  theme(text = element_text(size=12), 
        axis.text.x = element_text(size = 12, angle = 50, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

dev.off()

png(filename = "feedinghabit.png",width = 24, height = 20, units = "cm", pointsize = 12, res=600)

ggplot(gutdata, aes(x=feeding.habit , y=log(Mpsgut + 1))) + 
  geom_boxplot(size=1, fill = 'grey80') + 
  xlab("Feeding habit") + 
  ylab("ln(MPs + 1)") +
  theme_bw() +
  theme(text = element_text(size=12), 
        axis.text.x = element_text(size = 12, angle = 50, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

dev.off()

plot(N~feeding.habit, data = gutdata)
