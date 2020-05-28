library(plyr)
library(ggplot2)
library(nlme)
library(MuMIn)
library(lme4)
library(cowplot)
library(glmmTMB)
library(ggthemes)
library(colorspace)
library(DHARMa)
library(boot)

model.assess <- 
  function(x) {
  plot(resid(x, type = 'pearson') ~ fitted(x, type = 'response'))
    abline(0,0,lty=2)
  }

theme1 <-
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    legend.text = element_text(size = 7),
    panel.grid = element_blank()
  )

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
length(trophicfish2$species) # 873 data points
length(trophicfish$species) # consolidated from 977 data point
length(unique(trophicfish2$species)) # ~647 species
length(unique(trophicfish2$family)) # from 171 families
sum(na.omit(trophicfish2$N))  # 30,543 individuals  

trophicfish2$study <- with(trophicfish2, paste0(author, year))
trophicfish2$study <- as.factor(trophicfish2$study)
head(trophicfish2)
length(unique(trophicfish2$study))  # 124 + 1 studies

summary(trophicfish2$min.size)

hist(trophicfish2$min.size)

trophicfish2$polymer.ID <- 
  mapvalues(trophicfish2$polymer.ID, 
            from = c('yes','no'),
            to = c('Polymer ID Method Used',
                   'Polymer ID Method Not Used'))
trophicfish2$blanks <- 
  mapvalues(trophicfish2$blanks,
            from = c('yes', 'no'),
            to = c('Blanks Used',
                   'Blanks Not Used'))

trophicfish2$feeding.habit <- as.factor(trophicfish2$feeding.habit)

levels(trophicfish2$feeding.habit)

trophicfish2$feeding.habit <- mapvalues(trophicfish2$feeding.habit, 
                                   from = levels(trophicfish2$feeding.habit),
                                   to = c("Not listed", 
                                          "Browsing on substrate",
                                          "Detritus",
                                          "Filtering plankton",
                                          "Grazing on aquatic plants", 
                                          "Hunting macrofauna",
                                          "Other",
                                          "Selective plankton feeding", 
                                          "Variable"))

levels(trophicfish2$feeding.habit)

trophicfish2$feeding.habit <- factor(
  trophicfish2$feeding.habit,
  levels = c(
    "Not listed",
    "Detritus",
    "Browsing on substrate",
    "Grazing on aquatic plants",
    "Filtering plankton",
    "Selective plankton feeding",
    "Hunting macrofauna",
    "Variable",
    "Other"
  )
)

levels(trophicfish2$feeding.habit)

trophicfish2$study.habitat <- as.factor(trophicfish2$study.habitat)
trophicfish2$environment <- as.factor(trophicfish2$environment)

gutdata <- subset(trophicfish2, Mpsgut != 'NA')
summary(gutdata)
summary(gutdata$author)
length(gutdata$species) # 731 data points
length(unique(gutdata$species)) # 553 species
length(unique(gutdata$family)) # from 158 families
length(unique(gutdata$study)) # from 104 + 1 studies

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
  scale_y_continuous(trans = 'log1p') +
  theme_bw() +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(size = 16, angle = 50, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

hist(gutdata$Mpsgut)

## Trophic modeling

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

## Need to account for environment, region, and sample size

gutdata <- subset(gutdata, !is.na(TL) & !is.na(Mpsgut))  # remove NAs

gutdata$Mpsgut <- ceiling(gutdata$Mpsgut)  # round up to the nearest integer

mod1 <-
  glmmTMB(
    log(Mpsgut + 1) ~ TL + environment + (TL | region / study),
    weights = N,
    data = gutdata
  )

model.assess(mod1)  ## variance looks pretty good

plot(resid(mod1) ~ gutdata$TL)
plot(resid(mod1) ~ gutdata$environment)
plot(resid(mod1) ~ gutdata$region)
plot(resid(mod1) ~ gutdata$study)

summary(mod1)

mod1.1 <- glmmTMB(log(Mpsgut + 1) ~ TL + (TL | region / study), 
                  weights = N, data = gutdata)
anova(mod1, mod1.1)  # environment significant, p < 0.001


## Plot model predictions

ilink1 <- family(mod1)$linkinv
prediction1 <- predict(mod1, re.form = NA, se.fit = TRUE)

gutdata$predict <- exp(ilink1(prediction1$fit)) - 1
gutdata$upper <- exp(ilink1(prediction1$fit +
                              (1.96 * prediction1$se.fit))) - 1
gutdata$lower <- exp(ilink1(prediction1$fit -
                              (1.96 * prediction1$se.fit))) - 1

freshplot <-
  ggplot(subset(gutdata, study.habitat == 'Freshwater')) +
  geom_line(aes(x = TL, y = predict, colour = environment),
            size = 0.5, alpha = 0.8) +
  geom_ribbon(aes(x = TL, ymin = lower, ymax = upper, fill = environment), 
              alpha = 0.3) +
  geom_point(aes(x = TL, y = Mpsgut, colour = environment),
             shape = 1, size = 0.75) +
  facet_wrap(~ region, scales = 'free_x', ncol = 4,
             labeller = label_wrap_gen(width = 20)) +
  labs(x = 'Trophic Level',
       y = expression(paste(
         'Microplastic Concentration (particles ' ~ ind ^ -1 * ')'
       )),
       size = 'Sample Size') +
  coord_cartesian(xlim = c(2,5), ylim = c(0,30)) +
  scale_fill_manual(values = qualitative_hcl(n = 5, palette = 'Warm'),
                    name = 'Environment') +
  scale_colour_manual(values = qualitative_hcl(n = 5, palette = 'Warm'),
                      name = 'Environment') +
  scale_x_continuous(breaks = seq(from = 2, to = 5, by = 1)) +
  scale_y_continuous(trans = 'log1p', breaks = c(0, 1, 10, 30)) +
  theme1

marineplot <-
  ggplot(subset(gutdata, study.habitat == 'Marine')) +
  geom_line(aes(x = TL, y = predict, colour = environment),
            size = 0.5, alpha = 0.8) +
  geom_ribbon(aes(x = TL, ymin = lower, ymax = upper, fill = environment), 
              alpha = 0.3) +
  geom_point(aes(x = TL, y = Mpsgut, colour = environment),
             shape = 1, size = 0.75) +
  facet_wrap(~ region, scales = 'free_x', ncol = 4,
             labeller = label_wrap_gen(width = 20)) +
  labs(x = 'Trophic Level',
       y = expression(paste(
         'Microplastic Concentration (particles ' ~ ind ^ -1 * ')'
       )),
       size = 'Sample Size') +
  coord_cartesian(xlim = c(2,5), ylim = c(0,30)) +
  scale_fill_manual(values = qualitative_hcl(n = 8, palette = 'Cold'),
                    name = 'Environment') +
  scale_colour_manual(values = qualitative_hcl(n = 8, palette = 'Cold'),
                      name = 'Environment') +
  scale_x_continuous(breaks = seq(from = 2, to = 5, by = 1)) +
  scale_y_continuous(trans = 'log1p', breaks = c(0, 1, 10, 30)) +
  theme1


png('Gut Content Plot.png', width = 16.7, height = 20, units = 'cm', res = 300)

plot_grid(freshplot, marineplot, labels = c('A', 'B'), rel_heights = c(1,1.9),
          nrow = 2, align = 'v')

dev.off()


## Effect of region on ingestion rates
## Need to account for study FAO area

ingestion <- subset(trophicfish2, !is.na(IR))

summary(ingestion$IR)

length(ingestion$species) # 639 data points
length(unique(ingestion$species)) # 478 species
length(unique(ingestion$family)) # from 157 families
length(unique(ingestion$study)) # 105 studies

ingestion$region <- as.factor(ingestion$region)

levels(ingestion$region)

summary(ingestion$region)
ingestion$region <- as.character(ingestion$region)
ingestion$region <- as.factor(ingestion$region)

summary(ingestion)

ing.mod1 <- glmmTMB(
  IR ~ TL + (TL | region / study),
  weights = N,
  data = ingestion,
  family = beta_family(link = 'logit')
)

summary(ing.mod1)  ## TL not sig., p = 0.6
model.assess(ing.mod1) # variance looks pretty homogenous
plot(resid(ing.mod1, type = 'pearson') ~ 
       ingestion$region) # patterns, but not too wild

## Plot model predictions

ing.prediction <- predict(ing.mod1, type = 'response')
ing.upper <- ing.prediction + 
  (2*predict(ing.mod1, type = 'response', se.fit = TRUE)$se.fit)
ing.lower <- ing.prediction - 
  (2*predict(ing.mod1, type = 'response', se.fit = TRUE)$se.fit)

png('Ingestion Rate Plot.png', width = 33, height = 19, 
    units = 'cm', res = 300)

ggplot(ingestion) +
  geom_line(aes(x = TL, y = ing.prediction, colour = region),
            size = 1, alpha = 0.8) +
  geom_ribbon(aes(x = TL, ymin = ing.lower, ymax = ing.upper, fill = region, 
                  colour = region), 
              alpha = 0.3, linetype = 'dashed') +
  geom_jitter(aes(x = TL, y = IR, size = N, colour = region),
              shape = 1) +
  facet_grid(. ~ area, scales = 'free_x') +
  labs(x = 'Trophic Level',
       y = 'Ingestion Rate',
       colour = 'FAO Area',
       fill = 'FAO Area',
       size = 'Sample Size') +
  coord_cartesian(xlim = c(2,4.7), ylim = c(0,1)) +
  scale_x_continuous(breaks = seq(from = 2, to = 5, by = 0.5)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.1)) +
  theme_few() +
  scale_color_manual(values = col1) +
  scale_fill_manual(values = col1) +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 12)
  )

dev.off()

## Do microplastic numbers in the gut correlate with body size? 

allo <- subset(gutdata, total.length != 'NA')

length(allo$total.length)  # 342 data point remaining
length(unique(allo$study))  # 52 studies
length(unique(allo$species))  # 300 species

allo.mod1 <- 
  glmmTMB(log(Mpsgut + 1) ~ total.length + (total.length | region), 
          weights = N, data = allo)
summary(allo.mod1)
plot(resid(allo.mod1, type = 'pearson') ~ fitted(allo.mod1))  # not too bad
plot(resid(allo.mod1, type = 'pearson') ~ allo$total.length)
# variance appears to decrease with increasing total length

## try putting total length on a log scale

allo.mod2 <- 
  glmmTMB(log(Mpsgut + 1) ~ log(total.length) + (log(total.length) | region), 
          weights = N, data = allo)

model.assess(allo.mod2)  # not too bad
plot(resid(allo.mod2, type = 'pearson') ~ log(allo$total.length))  # better

AICc(allo.mod1, allo.mod2)  # allo.mod2 is a slightly better fit

summary(allo.mod2)  # total length not sig., p = 0.5

## Plot model predictions

allo.prediction <- exp(predict(allo.mod2, type = 'response')) - 1
allo.upper <- exp(predict(allo.mod2, type = 'response') + 
                    (2*predict(allo.mod2, type = 'response', 
                               se.fit = TRUE)$se.fit)) - 1
allo.lower <- exp(predict(allo.mod2, type = 'response') - 
                    (2*predict(allo.mod2, type = 'response', 
                               se.fit = TRUE)$se.fit)) - 1

png('Allometric Gut Concentration Plot.png', width = 33, height = 19, 
    units = 'cm', res = 300)

ggplot(allo) +
  geom_line(aes(x = log(total.length), y = allo.prediction, colour = region),
            size = 1, alpha = 0.8) +
  geom_ribbon(aes(x = log(total.length), ymin = allo.lower, ymax = allo.upper, fill = region, 
                  colour = region), 
              alpha = 0.3, linetype = 'dashed') +
  geom_point(aes(x = log(total.length), y = Mpsgut, size = N, colour = region),
              shape = 1) +
  facet_grid(. ~ area, scales = 'free_x') +
  labs(x = 'ln(Total Length) (cm)',
       y = expression(paste(
         'Microplastic Concentration (particles ' ~ ind ^ -1 * ')')),
       colour = 'FAO Area',
       fill = 'FAO Area',
       size = 'Sample Size') +
  coord_cartesian(ylim = c(0,30)) +
  scale_y_continuous(breaks = seq(from = 0, to = 30, by = 5),
                     expand = c(0,0), trans = 'log1p') +
  theme_few() +
  scale_color_manual(values = col1) +
  scale_fill_manual(values = col1) +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 12)
  )

dev.off()


## Does ingestion rate increase with body size?

allo2 <- subset(ingestion, total.length != 'NA')

allo2$species <- as.character(allo2$species)
allo2$species <- as.factor(allo2$species)

length(allo2$total.length)  # 247 data point remaining
length(unique(allo2$study))  #45 studies
length(unique(allo2$species))  #214 species

ingest.allo.mod1 <- 
  glmmTMB(IR ~ log(total.length) + (log(total.length) | region),
          family = binomial(), weights = N, data = allo2)
model.assess(ingest.allo.mod1)
summary(ingest.allo.mod1)  # p = 0.722

## Plot model predictions

ingest.allo.prediction <- predict(ingest.allo.mod1, type = 'response')
ingest.allo.upper <- ingest.allo.prediction + 
  (2*predict(ingest.allo.mod1, type = 'response', se.fit = TRUE)$se.fit)
ingest.allo.lower <- ingest.allo.prediction - 
  (2*predict(ingest.allo.mod1, type = 'response', se.fit = TRUE)$se.fit)

png('Allometric Ingestion Rate Plot.png', width = 33, height = 13, 
    units = 'cm', res = 300)

ggplot(allo2) +
  geom_line(aes(x = log(total.length), y = ingest.allo.prediction, 
                colour = region),
            size = 1, alpha = 0.8) +
  geom_ribbon(aes(x = log(total.length), 
                  ymin = ingest.allo.lower, 
                  ymax = ingest.allo.upper, 
                  fill = region, 
                  colour = region), 
              alpha = 0.3, linetype = 'dashed') +
  geom_jitter(aes(x = log(total.length), y = IR, size = N, colour = region),
              shape = 1) +
  facet_grid(. ~ area) +
  labs(x = 'ln(Total Length) (cm)',
       y = 'Ingestion Rate',
       colour = 'FAO Area',
       fill = 'FAO Area',
       size = 'Sample Size') +
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2),
                     expand = c(0,0)) +
  theme_few() +
  scale_color_manual(values = col1) +
  scale_fill_manual(values = col1) +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.box = 'vertical',
    legend.position = 'bottom'
  )

dev.off()


## Is there an effect of methodology on reported MP gut numbers?

methods <- subset(gutdata, 
                  min.size != 'NA' & polymer.ID != 'NA' &
                    blanks != 'NA')

meth.mod1 <- 
  glmmTMB(log(Mpsgut + 1) ~ scale(min.size, center = TRUE) + polymer.ID + 
            blanks + (scale(min.size, center = TRUE) | region), 
          weights = N, data = methods)

# try without varying intercept
meth.mod2 <- 
  glmmTMB(log(Mpsgut + 1) ~ scale(min.size, center = TRUE) + polymer.ID + 
            blanks + (1 | region), 
          weights = N, data = methods)

AICc(meth.mod1, meth.mod2) # better fit with varying intercept

model.assess(meth.mod2)

summary(meth.mod1)  # polymer ID and blanks sig., p < 0.001

## More MPs with polymer ID and with blanks
## But only on the order of ~0.1-0.3 particles

## Plot model predictions

meth.prediction <- exp(predict(meth.mod1, type = 'response')) - 1
meth.upper <- exp(predict(meth.mod1, type = 'response') + 
                    (2*predict(meth.mod1, type = 'response', 
                               se.fit = TRUE)$se.fit)) - 1
meth.lower <- exp(predict(meth.mod1, type = 'response') - 
                    (2*predict(meth.mod1, type = 'response', 
                               se.fit = TRUE)$se.fit)) - 1

png('Methods Plot.png', width = 27, height = 18, 
    units = 'cm', res = 300)

ggplot(methods) +
  geom_line(aes(x = min.size, y = meth.prediction, colour = region),
            size = 1, alpha = 0.8) +
  geom_ribbon(aes(x = min.size, ymin = meth.lower, ymax = meth.upper, 
                  fill = region, colour = region), 
              alpha = 0.3, linetype = 'dashed') +
  geom_point(aes(x = min.size, y = Mpsgut, size = N, colour = region),
             shape = 1) +
  facet_grid(polymer.ID ~ blanks) +
  labs(x = expression(paste('Lowest Detectible Particle Size ('*mu*'m)')),
       y = expression(paste('Microplastic Concentration (particles ' ~ 
                              ind ^ -1 * ')')),
       colour = 'FAO Area',
       fill = 'FAO Area',
       size = 'Sample Size') +
  coord_cartesian(ylim = c(0,40)) +
  scale_y_continuous(breaks = c(0,1,5,10,20,30,40),
                     expand = c(0,0), trans = 'log1p') +
  scale_x_continuous(breaks = c(0,1,10,100,500),
                     trans = 'log1p') +
  theme_few() +
  scale_color_manual(values = col1) +
  scale_fill_manual(values = col1) +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 12),
    panel.spacing = unit(c(0.5), units = 'cm')
  )

dev.off()

factors <- factor(c('Intercept', 'Lowest Detectable Particle Size', 
                    'Polymer ID Method Used', 'Blanks Used'))
coeff <- as.numeric(fixef(meth.mod1)[[1]])

trans.coeff <- exp(coeff) - 1

meth.coef <- data.frame(factors, coeff, trans.coeff)

meth.coef

plot(coeff ~ factors, data = meth.coef, ylim = c(-2,2))
abline(0,0, lty = 2)

## Try trophic level/regional model using max 100 micron LOD

gutdata100 <- subset(gutdata, min.size <= 100)

length(gutdata100$Mpsgut)  # 142 data points remaining

gutdata100$region <- as.character(gutdata100$region)
gutdata100$region <- as.factor(gutdata100$region)

gut2mod1 <- 
  glmmTMB(log(Mpsgut + 1) ~ TL + (TL | region), 
          weights = N, data = gutdata100)

model.assess(gut2mod1)

summary(gut2mod1)  # trophic level not sig., p = 0.9

## Plot model predictons

gut2.prediction <- exp(predict(gut2mod1))-1
gut2.upper <- exp(predict(gut2mod1) - 
               (2*predict(gut2mod1, se.fit = TRUE)$se.fit)) - 1
gut2.lower <- exp(predict(gut2mod1) + 
               (2*predict(gut2mod1, se.fit = TRUE)$se.fit)) - 1

png('Gut Content Reduced Plot.png', width = 33, height = 19, units = 'cm', 
    res = 300)

ggplot(subset(gutdata100, TL != 'NA')) +
  geom_line(aes(x = TL, y = gut2.prediction, colour = region),
            size = 1, alpha = 0.8) +
  geom_ribbon(aes(x = TL, ymin = gut2.lower, ymax = gut2.upper, fill = region, 
                  colour = region), 
              alpha = 0.3) +
  geom_point(aes(x = TL, y = Mpsgut, size = N, colour = region),
             shape = 1) +
  facet_grid(. ~ area, scales = 'free_x') +
  labs(x = 'Trophic Level',
       y = expression(paste(
         'Microplastic Concentration (particles ' ~ ind ^ -1 * ')'
       )),
       colour = 'FAO Area',
       fill = 'FAO Area',
       size = 'Sample Size') +
  coord_cartesian(xlim = c(2,4.7), ylim = c(0,30)) +
  scale_x_continuous(breaks = seq(from = 2, to = 5, by = 0.5)) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,10,15,20,25,30), 
                     trans = 'log1p', expand = c(0,0)) +
  theme_few() +
  scale_color_manual(values = col1) +
  scale_fill_manual(values = col1) +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 12)
  )

dev.off()



## Model in terms of concentration by body weight

gutdata$gutconc <- with(gutdata, Mpsgut/W)
summary(gutdata$gutconc)

gut.conc <- subset(gutdata, gutconc != 'NA' & TL != 'NA')
summary(gut.conc)

length(gut.conc$species) # 424 data points
length(unique(gut.conc$species)) # ~357 species
length(unique(gut.conc$family)) # from 110 families
length(unique(gut.conc$study)) # 55 studies
length(unique(gut.conc$region)) ## 18 regions

mag.mod1 <- 
  glmmTMB(gutconc ~ TL + (TL | region), weights = N, data = gut.conc)

model.assess(mag.mod1)

## Need to log transform gutconc

mag.mod2 <-
  glmmTMB(log(gutconc + 1) ~ TL + (TL | region), weights = N, data = gut.conc)

model.assess(mag.mod2)  ## Variance still looks odd, log-transform TL

mag.mod3 <-
  glmmTMB(log(gutconc + 1) ~ log(TL) + (log(TL) | region), weights = N, 
          data = gut.conc)

model.assess(mag.mod3)

AICc(mag.mod1, mag.mod2, mag.mod3)  # best as log-log

# try removing the random effect

mag.mod4 <- glmmTMB(log(gutconc + 1) ~ log(TL), weights = N, data = gut.conc)

model.assess(mag.mod4)

AICc(mag.mod3, mag.mod4)  # better fit with RE

## what's driving the variance pattern?

plot(resid(mag.mod2, type = 'pearson') ~ log(gut.conc$TL))
plot(resid(mag.mod2, type = 'pearson') ~ na.omit(gut.conc$region))

summary(mag.mod2)  # p = 0.158

## Plot model predictions

mag.prediction <- exp(predict(mag.mod2))-1
mag.lower <- exp(predict(mag.mod2) - 
                   (2*predict(mag.mod2, se.fit = TRUE)$se.fit)) - 1
mag.upper <- exp(predict(mag.mod2) + 
                   (2*predict(mag.mod2, se.fit = TRUE)$se.fit)) - 1

png('Bioaccumulation Plot.png', width = 33, height = 14, 
    units = 'cm', res = 300)

ggplot(gut.conc) +
  geom_line(aes(x = TL, y = mag.prediction, colour = region),
            size = 1, alpha = 0.8) +
  geom_ribbon(aes(x = TL, ymin = mag.lower, ymax = mag.upper, fill = region, 
                  colour = region), 
              alpha = 0.3, linetype = 'dashed') +
  geom_point(aes(x = TL, y = gutconc, size = N, colour = region),
             shape = 1) +
  facet_grid(. ~ area, scales = 'free_x') +
  labs(x = 'Trophic Level',
       y = expression(paste(
         'Microplastic Concentration (particles '*g^-1*')'
       )),
       colour = 'FAO Area',
       fill = 'FAO Area',
       size = 'Sample Size') +
  coord_cartesian(xlim = c(2,5), ylim = c(0,30)) +
  scale_x_continuous(breaks = seq(from = 2, to = 5, by = 1)) +
  scale_y_continuous(breaks = c(0,1,5,10,20,30), 
                     trans = 'log1p', expand = c(0,0)) +
  theme_few() +
  scale_color_manual(values = col1) +
  scale_fill_manual(values = col1) +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.box = 'vertical',
    legend.position = 'bottom'
  )

dev.off()

png('Bioaccumulation Plot Combined.png', width = 15, height = 13, 
    units = 'cm', res = 300)

ggplot(gut.conc) +
  geom_line(aes(x = TL, y = mag.prediction, group = region),
            size = 1, alpha = 0.8) +
  geom_point(aes(x = TL, y = gutconc, size = N),
             shape = 1) +
  labs(x = 'Trophic Level',
       y = expression(paste(
         'Microplastic Concentration (particles '*g^-1*')'
       )), size = 'Sample Size') +
  coord_cartesian(xlim = c(2,5), ylim = c(0,30)) +
  scale_x_continuous(breaks = seq(from = 2, to = 5, by = 1)) +
  scale_y_continuous(breaks = c(0,1,5,10,20,30), 
                     trans = 'log1p', expand = c(0,0)) +
  theme_few() +
  scale_color_manual(values = col1) +
  scale_fill_manual(values = col1) +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.box = 'vertical',
    legend.position = 'bottom'
  )

dev.off()



## Plot according to lower limit of detection

png('Lower Limit Plot.png', width = 30, height = 19, 
    units = 'cm', res = 300)

ggplot(gutdata2) +
  geom_point(aes(
    x = reorder(region, Mpsgut, mean),
    y = Mpsgut,
    size = N,
    colour = min.size
  ),
  alpha = 0.5) +
  labs(
    x = 'Region',
    y = expression(paste(
      'Microplastic Concentration (particles ' ~
        ind ^ -1 * ')'
    )),
    colour = expression(paste('Lowest Detectable Particle Size ('*mu*'m)')),
    size = 'Sample Size'
  ) +
  coord_cartesian(ylim = c(0, 35)) +
  scale_y_continuous(
    breaks = c(0, 1, 5, 10, 20, 30),
    expand = c(0, 0),
    trans = 'log1p'
  ) +
  coord_flip() +
  scale_colour_continuous_sequential(palette = 'Red-Blue') +
  theme_few() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 12),
    panel.spacing = unit(c(0.5), units = 'cm')
  )

dev.off()
