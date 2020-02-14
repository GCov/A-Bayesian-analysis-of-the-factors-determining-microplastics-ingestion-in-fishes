library(plyr)
library(ggplot2)
library(nlme)
library(MuMIn)
library(lme4)
library(cowplot)
library(car)
library(mgcv)
library(glmmTMB)
library(ggthemes)
library(merTools)

trophicfish <- 
  read.csv("trophicfishdata.csv", header = TRUE)

names(trophicfish)
trophicfish <- trophicfish[1:42]

summary(trophicfish)

trophicfish$year <- as.factor(trophicfish$year)

# concentration by weight

trophicfish$conc <- trophicfish$Mpsgut/trophicfish$W
summary(trophicfish)
trophicfish$IR <- trophicfish$ingest.rate/100

## Consolidate averages by study

trophicfish$fork.lengthprod <- trophicfish$fork.length*trophicfish$N
trophicfish$total.lengthprod <- trophicfish$total.length*trophicfish$N
trophicfish$Wprod <- trophicfish$W*trophicfish$N
trophicfish$GWprod <- trophicfish$GW*trophicfish$N
trophicfish$IRprod <- trophicfish$IR*trophicfish$N
trophicfish$Mpsgutprod <- trophicfish$Mpsgut*trophicfish$N
trophicfish$ SDMPsgutprod <- trophicfish$ SDMPsgut*trophicfish$N

trophicfish2 <- ddply(trophicfish, 
                      c('author', 'study.habitat', 'year', 'region', 'species', 
                        'family', 'genus', 'environment', 'climate', 'red.list', 
                        'feeding.type', 'feeding.habit', 'TL', 'min.size', 
                        'float.meth', 'dig.meth', 'count.meth', 'polymer.meth',
                        'polymer.ID', 'blanks', 'maj.fib', 'maj.under.one.mm', 
                        'maj.polymer', 'maj.col', 'exclude.fib'), 
                      summarise, 
                      fork.lengthprod = sum(fork.lengthprod), 
                      total.lengthprod = sum(total.lengthprod), 
                      Wprod = sum(Wprod), 
                      GWprod = sum(GWprod),
                      N = sum(N), 
                      IRprod = sum(IRprod), 
                      Mpsgutprod = sum(Mpsgutprod), 
                      SDMPsgutprod = sum(SDMPsgutprod)
)

trophicfish2$fork.length <- trophicfish2$fork.lengthprod/trophicfish2$N
trophicfish2$total.length <- trophicfish2$total.lengthprod/trophicfish2$N
trophicfish2$W <- trophicfish2$Wprod/trophicfish2$N
trophicfish2$GW <- trophicfish2$GWprod/trophicfish2$N
trophicfish2$IR <- trophicfish2$IRprod/trophicfish2$N
trophicfish2$Mpsgut <- trophicfish2$Mpsgutprod/trophicfish2$N
trophicfish2$SDMPsgut <- trophicfish2$SDMPsgutprod/trophicfish2$N


summary(trophicfish2)
length(trophicfish2$species) # 780 dta points
length(trophicfish$species) # consolidated from 872 data point
length(unique(trophicfish2$species)) # ~607 species
length(unique(trophicfish2$family)) # from 165 families

trophicfish2$study <- with(trophicfish2, paste0(author, year))
trophicfish2$study <- as.factor(trophicfish2$study)
head(trophicfish2)
length(unique(trophicfish2$study))  # 109 + 1 studies

summary(trophicfish2$min.size)

hist(trophicfish2$min.size)

gutdata <- subset(trophicfish2, Mpsgut != 'NA')
summary(gutdata)
summary(gutdata$author)
length(gutdata$species) # 654 data points
length(unique(gutdata$species)) # 521 species
length(unique(gutdata$family)) # from 152 families
length(unique(gutdata$study)) # from 91 + 1 studies

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

ggplot(gutdata, aes(x=climate , y=Mpsgut)) + 
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

gutdata$gutconc <- with(gutdata, Mpsgut/W)
summary(gutdata$gutconc)

gut.conc <- subset(gutdata, gutconc != 'NA')
summary(gut.conc)

length(gut.conc$species) # 428 data points
length(unique(gut.conc$species)) # ~361 species
length(unique(gut.conc$family)) # from 111 families
length(unique(gut.conc$study)) # 55 studies
length(unique(gut.conc$region)) ## 20 regions

levels(gutdata$feeding.habit)

gutdata$feeding.habit <- mapvalues(gutdata$feeding.habit, 
                                   from = levels(gutdata$feeding.habit),
                                   to = c("Not listed", 
                                          "Browsing on substrate", 
                                          "Filtering plankton",
                                          "Grazing on aquatic plants", 
                                          "Hunting macrofauna",
                                          "Other",
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

## Does # MPs in the gut correlate with trophic level?

## Need to account for freshwater vs. marine, region, and sample size

mod1 <- glmmTMB(Mpsgut ~ TL*study.habitat + (TL | region), 
                weights = N, data = gutdata)
summary(mod1)
plot(resid(mod1, type = 'pearson') ~ fitted(mod1))  ## variance is a bit weird

# try log transforming Mpsgut

mod2 <- glmmTMB(log(Mpsgut + 1) ~ TL*study.habitat + (TL | region), 
                weights = N, data = gutdata)
summary(mod2)
plot(resid(mod2, type = 'pearson') ~ fitted(mod2)) 
plot(resid(mod2, type = 'pearson') ~ na.omit(gutdata$TL))
plot(resid(mod2, type = 'pearson') ~ subset(gutdata, TL != 'NA')$study.habitat)
AICc(mod1, mod2) # way better model fit, variance doesn't look too bad

## Plot model predictions

prediction <- exp(predict(mod2)-1)
upper <- exp(predict(mod2) - 
               (2*predict(mod2, se.fit = TRUE)$se.fit) - 1)
lower <- exp(predict(mod2) + 
               (2*predict(mod2, se.fit = TRUE)$se.fit) - 1)

png('Gut Content Plot.png', width = 30, height = 15, units = 'cm', res = 300)

ggplot(subset(gutdata, TL != 'NA')) +
  geom_line(aes(x = TL, y = prediction, colour = region),
            size = 1.5) +
  geom_ribbon(aes(x = TL, ymin = lower, ymax = upper, fill = region, 
                  colour = region), 
              alpha = 0.5) +
  geom_jitter(aes(x = TL, y = Mpsgut, size = N, colour = region),
              shape = 1) +
  facet_grid(. ~ study.habitat) +
  labs(x = 'Trophic Level',
       y = expression(paste(
         'Microplastic Concentration (particles ' ~ ind ^ -1 * ')'
       ))) +
  coord_cartesian(xlim = c(2,4.7), ylim = c(0,30)) +
  scale_x_continuous(breaks = seq(from = 2, to = 5, by = 0.5)) +
  scale_y_continuous(breaks = seq(from = 0, to = 30, by = 5), 
                     trans = 'log1p') +
  theme_few() +
  scale_color_viridis_d(option = 'C') +
  scale_fill_viridis_d(option = 'C') +
  guides(color = FALSE, fill = FALSE) +
  
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 16)
  )

dev.off()

ggplot(gutdata, aes(x = reorder(region, predict(mod2, level = 0), mean))) + 
  geom_errorbar(aes(ymin = predict(mod2, level = 0) - 
                    predict(mod2, level = 0, se.fit = TRUE)$se.fit,
                  ymax = predict(mod2, level = 0) + 
                    predict(mod2, level = 0, se.fit = TRUE)$se.fit
                  )) +
  geom_point(aes(y = predict(mod2, level = 0))) +
  labs(y = '# MPs in Gut', x = 'FAO Region', size = 'Number of Studies') +
  coord_flip() +
  theme_classic() + 
  theme(text = element_text(size=7), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7))

coef(mod4)

## Effect of region on ingestion rates
## Need to account for study FAO area

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

summary(ingestion)

mod7 <- glm(IR ~ TL + region + study,
              weights = N,
              data = ingestion, 
              family = binomial(link = 'logit'))
summary(mod7)  ## Model will not converge with study as a random effect
plot(resid(mod7) ~ fitted(mod7))

drop1(mod7)

ggplot(ingestion) +
  geom_point(aes(x = TL,
                 y = IR),
             size = 0.5) +
  geom_line(aes(x = TL,
                 y = predict(mod7, type = 'response', level = 0),
                 colour = study)) +
  labs(y = 'Ingestion Rate', x = 'Trophic Level') +
  facet_wrap(~ region) +
  theme_classic() + 
  theme(text = element_text(size=7), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = 'none')

## Do microplastic numbers in the gut correlate with body size? 

allo <- subset(gutdata, total.length != 'NA' & W != 'NA')

length(allo$total.length)  # 146 data point remaining

allo$cen.W <- (mean(allo$W) - allo$W)/sd(allo$W)
allo$cen.tot.l <- (mean(allo$total.length) - allo$total.length)/sd(allo$total.length)


mod8 <- lme(
  log(Mpsgut + 1) ~ cen.W + cen.tot.l,
  random = ~ 1 | study,
  weights = varPower(form =  ~ N),
  data = allo,
  method = "REML"
)

summary(mod8)
plot(resid(mod8) ~ fitted(mod8))  # variance does not look homogenous

## try without variance structure

mod9 <- lme(
  log(Mpsgut + 1) ~ cen.W + cen.tot.l,
  random = ~ 1 | study,
  data = allo,
  method = "REML"
)

plot(resid(mod9) ~ fitted(mod9))  # doesn't look much better

AICc(mod8, mod9)  # fit is a little better withvariance structure

# figure out what is driving the increasing variance

plot(resid(mod8) ~ allo$cen.W)
plot(resid(mod8) ~ allo$cen.tot.l)  # looks like it's total length and weight

# try log-transformaing total length and total width before centering and scaling

allo$cen.log.W <- (mean(log(allo$W)) - log(allo$W)) / sd(log(allo$W))
allo$cen.log.tot.l <- (mean(log(allo$total.length)) -
                         log(allo$total.length)) /
  sd(log(allo$total.length))


mod10 <- lme(log(Mpsgut+1) ~ cen.log.W + cen.log.tot.l, random = ~1|study, 
          weights = varPower(form =~ N), 
          data = allo, method = "REML")

plot(resid(mod10) ~ fitted(mod10))  # looks better

AICc(mod8, mod10)  # slightly better AICc value

plot(resid(mod10) ~ allo$cen.log.tot.l)
plot(resid(mod10) ~ allo$cen.log.W)  # way better

drop1(mod10)

mod11 <- lme(log(Mpsgut+1) ~ cen.log.W + cen.log.tot.l + region, random = ~1|study, 
             weights = varPower(form =~ N), 
             data = allo, method = "ML")

drop1(mod11)

mod12 <- update(mod11, ~ . -cen.log.tot.l)

AICc(mod11, mod12)

drop1(mod12)  # can remove the effect of region

mod13 <- update(mod12, ~. -region)

AICc(mod12, mod13)

summary(mod13)

plot(exp(predict(mod13, type = 'response')) ~ log(allo$W))

## number of MPs in gut appear to increase a bit with increasing weight
## but not a whole lot


## Does ingestion rate increase with body size?

allo2 <- subset(ingestion, total.length != 'NA' & W != 'NA')

length(allo2$total.length)  # 130 data point remaining

allo2$cen.log.W <- (mean(log(allo2$W)) - log(allo2$W)) / sd(log(allo2$W))
allo2$cen.log.tot.l <- (mean(log(allo2$total.length)) -
                         log(allo2$total.length)) /
  sd(log(allo2$total.length))

mod14 <- glm(IR ~ cen.log.tot.l + cen.log.W + region,
            weights = N,
            data = allo2, 
            family = binomial(link = 'logit'))
summary(mod14)

plot(resid(mod14) ~ fitted(mod14))

drop1(mod14)  ## best fitting model, can't remove anything

plot(predict(mod14, type = 'response') ~ allo2$cen.log.tot.l)
plot(predict(mod14, type = 'response') ~ allo2$cen.log.W)
plot(predict(mod14, type = 'response') ~ allo2$region)


## Is there an effect of methodology on reported MP gut numbers?

mod15 <- lmer(log(Mpsgut + 1) ~ min.size + polymer.ID + blanks + 
                region + (1 | study),
              data = gutdata)

plot(resid(mod15) ~ fitted(mod15))

drop1(mod15)  # remove region

mod16 <- update(mod15, ~ . -region)

AICc(mod15, mod16)

drop1(mod16)  # can remove polymer.ID

mod17 <- update(mod16, ~ . -polymer.ID)
AICc(mod16, mod17)

drop1(mod17)  # can remove blanks

mod18 <- update(mod17, ~ . -blanks)

AICc(mod17, mod18)

drop1(mod18)  # stop here

summary(mod18)

plot(resid(mod18) ~ fitted(mod18))

plot(resid(mod18) ~ gutdata$min.size)

plot(exp(predict(mod18, type = 'response'))-1 ~ gutdata$min.size)

plot(min.size ~ region, data = gutdata)

png(
  filename = "LOD Plot.png",
  width = 17,
  height = 17,
  units = "cm",
  pointsize = 7,
  res = 600
)

ggplot(gutdata) +
  geom_jitter(aes(x = reorder(region, Mpsgut, mean),
                 y = Mpsgut,
                 fill = min.size),
             colour = 'black', size = 1.5, shape = 21, alpha = 0.5) +
  labs(y = '# MPs in Gut', x = 'FAO Region', fill = 'Limit of Detection') +
  scale_y_continuous(trans = 'log1p') +
  coord_flip() +
  theme_classic() + 
  scale_colour_distiller(type = 'div', 
                         palette = 'RdYlBu', 
                         aesthetics = 'fill') +
  theme(text = element_text(size=7), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7))

dev.off()

## Try trophic level/regional model using max 100 micron LOD

gutdata100 <- subset(gutdata, min.size <= 100)

length(gutdata100$Mpsgut)  # 142 data points remaining

gutdata100$region <- as.character(gutdata100$region)
gutdata100$region <- as.factor(gutdata100$region)

plot(Mpsgut ~ region, data = gutdata100)
plot(Mpsgut ~ TL, data = gutdata100)

mod19 <- lme(
  log(Mpsgut + 1) ~ TL + environment + climate + region,
  weights = varExp(form = ~ N),
  random = ~ 1 | study,
  method = 'ML',
  data = gutdata100
)

plot(resid(mod19, type = 'pearson') ~ predict(mod19, type = 'response'))
plot(resid(mod19, type = 'pearson') ~ gutdata100$N)

## Try without variance structure

mod20 <- lme(
  log(Mpsgut + 1) ~ TL + environment + climate + region,
  random = ~ 1 | study,
  method = 'ML',
  data = gutdata100
)

AICc(mod19, mod20)  # better fit without variance structure

drop1(mod20)  # remove environment

mod21 <- update(mod20, ~. -environment)

drop1(mod21)  # remove climate

mod22 <- update(mod21, ~. -climate)

drop1(mod22)  # stop here
summary(mod22)

plot(resid(mod22, type = 'pearson') ~ fitted(mod22, type = 'response'))
plot(resid(mod22) ~ gutdata100$region)  ## pattern in residuals by region
plot(resid(mod22) ~ gutdata100$TL)

mod23 <- update(mod22, ~. -region)
mod24 <- update(mod22, ~. -TL)

anova(mod22, mod23)  # region significant, p = 0.01
anova(mod22, mod24)  # trophic level significant, p = 0.02

ggplot(gutdata100) +
  geom_ribbon(aes(x = TL,
                  ymin = exp(predict(mod22, type = 'response', level = 0) - 
                               predict(mod22, type = 'response', level = 0,
                                       se.fit = TRUE)$se.fit),
                  ymax = exp(predict(mod22, type = 'response', level = 0) + 
                               predict(mod22, type = 'response', level = 0,
                                       se.fit = TRUE)$se.fit),
                  fill = region
                  ),
              alpha = 0.5) +
  geom_line(aes(x = TL,
                y = exp(predict(
                  mod22, type = 'response', level = 0)),
                colour = region
                ),
            linetype = 'dashed') +
  geom_jitter(aes(x = TL,
                 y = Mpsgut,
                 colour = region),
             size = 0.5) +
  labs(y = 'Ingestion Rate', x = 'Trophic Level') +
  scale_y_continuous(trans = 'log1p') +
  scale_color_brewer(type = 'qual', palette = 'Set3', 
                     aesthetics = c('fill', 'colour')) +
  theme_classic() +
  theme(
    text = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    legend.position = 'none'
  )




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

ggplot(gut.conc, aes(x=TL , y=gutconc)) +
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






#################3

M1 <- lme(log(Mpsgut + 1) ~ TL, random = ~1|author/year, data = gutdata, 
          method = "REML")
    summary(M1)
plot(resid(M1) ~ fitted(M1))  # weird pattern for lower concentrations
plot(resid(M1) ~ gutdata$TL)  # looks good
plot(resid(M1) ~ gutdata$N)  # large variation according to sample size

abline(0,0)

## Compare fit with variance structure

M2 <- lme(log(Mpsgut+1) ~ TL, random = ~1|study, 
          weights = var(form =~ N), 
          data = gutdata, method = "REML")

AICc(M1, M2)  # better fit without variance structure

plot(resid(M2) ~ fitted(M2))

plot(resid(M1) ~ gutdata$TL)
plot(resid(M1) ~ gutdata$N)

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


M18 <- glm(IR ~ region + environment + TL, 
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
