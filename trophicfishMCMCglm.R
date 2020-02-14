library(plyr)
library(ggplot2)
library(nlme)
library(MuMIn)
library(lme4)
library(cowplot)
library(MCMCglmm)

trophicfish <- 
  read.csv("~/GRAD school/Trophic Fish Paper/Data/trophicfishupdate.csv")

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
                      c('author', 'year', 'region', 'species', 'family', 
                        'genus', 'environment', 'climate', 'red.list', 
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
length(trophicfish2$species) # 480 dta points
length(trophicfish$species) # consolidated from 557 data point
length(unique(trophicfish2$species)) # 383 species
length(unique(trophicfish2$family)) # from 131 families

trophicfish2$study <- with(trophicfish2, paste0(author, year))
trophicfish2$study <- as.factor(trophicfish2$study)
head(trophicfish2)
length(unique(trophicfish2$study))  # 63 studies

summary(trophicfish2$min.size)

summary(trophicfish2$min.size)

hist(trophicfish2$min.size)

trophicfish2$climate <- 
  mapvalues(trophicfish2$climate, from = "Benthopelagic", to = "Temperate")

gutdata <- subset(trophicfish2, Mpsgut != 'NA')
summary(gutdata)
summary(gutdata$author)
length(gutdata$species) # 366 data points
length(unique(gutdata$species)) # 384 species
length(unique(gutdata$family)) # from 121 families
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

length(gut.conc$species) # 194 data points
length(unique(gut.conc$species)) # 171 species
length(unique(gut.conc$family)) # from 70 families
length(unique(gut.conc$study)) # 25 studies
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

plot(Mpsgut ~ N, data = gutdata)  # variance decreases with N
plot(N ~ TL, data = gutdata)

gutdata$negN <- 0 - gutdata$N

plot(log(Mpsgut + 1) ~ negN, data = gutdata)  # seems to solve the issue

plot(Mpsgut ~ TL, data = gutdata) # variance doesn't look too bad

### MCMCglmm glm stuff

MCMCmod1 <- MCMCglmm(Mpsgut ~ TL + environment + climate + region, 
                random = ~ study,
                rcov = ~ N,
                family = NULL,
                data = gutdata)
