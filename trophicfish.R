summary(trophicfish)

plot(Mpsgut ~ TL, data = trophicfish)
plot(ingest.rate ~ TL, data = trophicfish)

plot(Mpsgut ~ region, data = trophicfish)
plot(ingest.rate ~ region, data = trophicfish)

plot(Mpsgut ~ family, data = trophicfish)
plot(ingest.rate ~ family, data = trophicfish)

plot(Mpsgut ~ environment, data = trophicfish)
plot(ingest.rate ~ environment, data = trophicfish)

plot(Mpsgut ~ climate, data = trophicfish)
plot(ingest.rate ~ climate, data = trophicfish)

plot(Mpsgut ~ red.list, data = trophicfish)
plot(ingest.rate ~ red.list, data = trophicfish)

plot(Mpsgut ~ feeding.type, data = trophicfish)
plot(ingest.rate ~ feeding.type, data = trophicfish)

plot(Mpsgut ~ feeding.habit, data = trophicfish)
plot(ingest.rate ~ feeding.habit, data = trophicfish)

plot(Mpsgut ~ TL.1, data = trophicfish)
plot(ingest.rate ~ TL.1, data = trophicfish)

plot(Mpsgut ~ W, data = trophicfish)
plot(ingest.rate ~ W, data = trophicfish)

plot(Mpsgut ~ min.size, data = trophicfish)
plot(ingest.rate ~ min.size, data = trophicfish)

plot(Mpsgut ~ polymer.meth, data = trophicfish)
plot(ingest.rate ~ polymer.meth, data = trophicfish)

head(trophicfish)

library(dplyr)


# concentration by weight

trophicfish$conc <- trophicfish$Mpsgut/trophicfish$W
plot(conc~TL, data = trophicfish) 
plot(Mpsgut~TL, data = trophicfish)
plot(conc~feeding.habit, data = trophicfish)
plot(conc~region, data = trophicfish)
plot(conc~environment, data = trophicfish)


summary(trophicfish)
trophicfish$IR <- trophicfish$ingest.rate/100

library(plyr)

## Consolidate averages by study

trophicfish$FLprod <- trophicfish$FL*trophicfish$N
trophicfish$TL.1prod <- trophicfish$TL.1*trophicfish$N
trophicfish$Wprod <- trophicfish$W*trophicfish$N
trophicfish$GWprod <- trophicfish$GW*trophicfish$N
trophicfish$IRprod <- trophicfish$IR*trophicfish$N
trophicfish$Mpsgutprod <- trophicfish$Mpsgut*trophicfish$N
trophicfish$ SDMPsgutprod <- trophicfish$ SDMPsgut*trophicfish$N

trophicfish2 <- ddply(trophicfish, c('author', 'year', 'region', 'species', 'family', 'genus', 'environment',
                                     'climate', 'feeding.type', 'feeding.habit', 'TL', 'min.size', 'float.meth',
                                     'dig.meth', 'count.meth', 'polymer.meth', 'maj.fib', 'maj.under.one.mm'), 
                      summarise, 
                      FLprod=sum(FLprod), TL.1prod=sum(TL.1prod), Wprod=sum(Wprod), GWprod=sum(GWprod),
                      N=sum(N), IRprod=sum(IRprod), Mpsgutprod=sum(Mpsgutprod), 
                      SDMPsgutprod=sum(SDMPsgutprod)
                      )

trophicfish2$FL <- trophicfish2$FLprod/trophicfish2$N
trophicfish2$TL.1 <- trophicfish2$TL.1prod/trophicfish2$N
trophicfish2$W <- trophicfish2$Wprod/trophicfish2$N
trophicfish2$GW <- trophicfish2$GWprod/trophicfish2$N
trophicfish2$IR <- trophicfish2$IRprod/trophicfish2$N
trophicfish2$Mpsgut <- trophicfish2$Mpsgutprod/trophicfish2$N
trophicfish2$ SDMPsgut <- trophicfish2$SDMPsgutprod/trophicfish2$N

summary(trophicfish2)
length(trophicfish2$species) # 370 dta points
length(trophicfish$species) # consolidated from 415 data point
length(unique(trophicfish2$species)) # 308 species
length(unique(trophicfish2$family)) # from 114 families
length(unique(trophicfish2$author)) # from 40 studies



summary(trophicfish2$min.size)

library(plyr)

trophicfish2$min.size <- mapvalues(trophicfish2$min.size, from = c("", "None"), 
                                   to = c(0, 0))

summary(trophicfish2$min.size)

trophicfish2$min.size <- as.numeric(levels(trophicfish2$min.size))[trophicfish2$min.size]

summary(trophicfish2$min.size)
hist(trophicfish2$min.size)


plot(Mpsgut ~ TL, data = trophicfish2)
plot(Mpsgut ~ region, data = trophicfish2)
plot(Mpsgut ~ environment, data = trophicfish2)

trophicfish2$climate <- mapvalues(trophicfish2$climate, from = "Benthopelagic", to = "Temperate")
  
library(MuMIn)

gutdata <- subset(trophicfish2, Mpsgut != 'NA')
summary(gutdata)
summary(gutdata$author)
length(gutdata$species) # 269 data points
length(unique(gutdata$species)) # 235 species
length(unique(gutdata$family)) # from 105 families
length(unique(gutdata$author)) # from 32 studies
tapply(gutdata$year, gutdata$author, unique)


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

length(gut.conc$species) # 131 data points
length(unique(gut.conc$species)) # 123 species
length(unique(gut.conc$family)) # from 60 families
length(unique(gut.conc$author)) # 16 studies
length(unique(gut.conc$region)) ## 9 regions
tapply(gut.conc$year, gut.conc$author, unique)

summary(gut.conc)


library(nlme)


gutdata$feeding.habit <- mapvalues(gutdata$feeding.habit, from = levels(gutdata$feeding.habit),
                                   to = c("Not listed", "Browsing on substrate", "Filtering plankton",
                                          "Grazing on aquatic plants", "Hunting macrofauna",
                                          "Selective plankton feeding", "Variable"))

levels(gutdata$feeding.habit)

gutdata$feeding.habit <- factor(gutdata$feeding.habit, levels = c("Not listed", "Browsing on substrate", 
                                                                  "Grazing on aquatic plants",
                                                                  "Filtering plankton", 
                                                                  "Selective plankton feeding", 
                                                                  "Hunting macrofauna",
                                                                  "Variable"))
levels(gutdata$feeding.habit)



hist(gutdata$Mpsgut)

plot(Mpsgut ~ N, data = gutdata) # variance definitely decreases with increasing N

M1 <- lme(log(Mpsgut+1) ~ TL + feeding.habit, random = ~1 + 1|author/year, weights = varExp(form =~ N), 
          data = gutdata, method = "ML")
summary(M1)
plot(M1)
AICc(M1)
gutdata$feeding.habit <- relevel(gutdata$feeding.habit, "Hunting macrofauna")

plot(resid(M1) ~ gutdata$feeding.habit)
plot(resid(M1) ~ log(gutdata$TL))

M1.1 <- lme(log(Mpsgut+1) ~ TL, random = ~1 + 1|author/year, weights = varExp(form =~ N), 
          data = gutdata, method = "ML")
AICc(M1, M1.1)
anova(M1, M1.1) # feeding habit is not significant, p<0.02

M1.2 <- lme(log(Mpsgut+1) ~ feeding.habit, random = ~1 + 1|author/year, weights = varExp(form =~ N), 
          data = gutdata, method = "ML")
AICc(M1, M1.2)
anova(M1, M1.2) # trophic level is significant, p<0.01


plot(log(Mpsgut+1) ~ log(TL), data = gutdata)
lines(fitted(M1)[sort(gutdata$TL)] ~ log(sort(gutdata$TL)), col = "red")

## Make predictions

TL <- seq(from = 2, to = 5, by = 0.1)
feeding.habit <- rep("Hunting macrofauna", length(TL))
author <- rep(0, length(TL))
N <- rep(100, length(TL))
predict.M1 <- data.frame(TL,feeding.habit, author, N)
head(predict.M1)


predict.M1$prediction <- predict(M1, predict.M1, level = 0)
head(predict.M1)

plot(prediction ~ TL, data = predict.M1)



## Account for weight

BSM.2 <- lm(log(W) ~ TL,  data = gut.conc)
AICc(BSM.2)
plot(BSM.2)
summary(BSM.2) ## log(W) loosely related to trophic level by log(y) = 1.25(x) -0.28

summary(gutdata$feeding.habit)

coefTable(M1)

gutdata$M1predict <- fitted(M1)

av.MPs <- gutdata
av.MPs$author <- 0
av.MPs$year <- 0
head(av.MPs)
gutdata$predict.av <- predict(M1, av.MPs, level = 0)
head(gutdata)

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


####################################
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
