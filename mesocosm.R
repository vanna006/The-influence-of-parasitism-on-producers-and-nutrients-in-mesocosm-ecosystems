### R code associated with Vannatta and Minchella. Food Webs 2021

###set directory

setwd("E:/A_in_arev/Microcosms/Data")


###import data
mydata <- read.csv("mesocosm.master.csv", header = TRUE)


### mean infection intensity in species
phymet <- lm(log(mydata$physa.mean.metas) ~ mydata$treatment)
summary(phymet)
plot(phymet)
shapiro.test(phymet$residuals)

promet <- lm(log(mydata$pro.mean.metas + 1) ~ mydata$treatment)
summary(promet)
plot(promet)
shapiro.test(promet$residuals)



###prevalence functions
phyprev <- lm(asinh(mydata$physa.prev) ~ mydata$treatment)
summary(phyprev)
shapiro.test(phyprev$residuals)

library(bestNormalize)
bestNormalize(phyprev$residuals)

library(mblm)
model.k = mblm(physa.prev ~ treatment,
               data=mydata)
summary(model.k)

proprev <- lm(mydata$pro.prev ~ mydata$treatment)
summary(proprev)
shapiro.test(proprev$residuals)

#########################
####### Analyses using mean infection intensity in Physa as a predictor
#########################
physa <- lm(log(mydata$physa.abun) ~ mydata$physa.mean.metas)
summary(physa)
shapiro.test(physa$residuals)#significance

promen <- lm(mydata$pro.abun ~ mydata$physa.mean.metas)
summary(promen)
shapiro.test(promen$residuals)

total.ab <- lm(log(mydata$total.abun) ~ mydata$physa.mean.metas)
summary(total.ab)
shapiro.test(total.ab$residuals)#significance

phys.pro <- lm(log(mydata$phys.pro) ~ mydata$physa.mean.metas)
summary(phys.pro)
shapiro.test(phys.pro$residuals)

phys.all <- lm((mydata$physa.abun/mydata$total.abun) ~ mydata$physa.mean.metas)
summary(phys.all)
shapiro.test(phys.all$residuals)


### periphyton characteristics
#create quadratic value for infected-snail days
infected2 <- mydata$physa.mean.metas * mydata$physa.mean.metas

#periphyton percent ash-free dry mass
AFDMq <-lm(mydata$X.AFDM ~ mydata$physa.mean.metas)
summary(AFDMq)
plot(AFDMq)
shapiro.test(AFDMq$residuals) #significance

### chlorophyll a
chloro <- lm(mydata$chla.per.cm2 ~ mydata$physa.mean.metas)
summary(chloro)
plot(chloro)
shapiro.test(chloro$residuals)

#periphyton dry mass per square centimeter
periq <-lm(log(mydata$dm.per.cm2) ~ mydata$physa.mean.metas)# + infected2)
summary(periq)
shapiro.test(periq$residuals)
plot(periq)#significance



###Wolfiella biomass
wolf.mod <- lm(sqrt(mydata$total.wolf) ~ mydata$physa.mean.metas)
summary(wolf.mod)
plot(wolf.mod)
shapiro.test(wolf.mod$residuals)

surfveg.mod <- lm(mydata$total.mass ~ mydata$physa.mean.metas)
summary(surfveg.mod)
shapiro.test(surfveg.mod$residuals)

lemna.mod <- lm(log(mydata$total.lemna) ~ mydata$physa.mean.metas)
summary(lemna.mod)
shapiro.test(lemna.mod$residuals)

wolf.lem <- lm(log(mydata$wolf.lem) ~ mydata$physa.mean.metas)
summary(wolf.lem)
shapiro.test(wolf.lem$residuals)

###in situ diel primary production
#week 6
DO <- lm(log(mydata$do.delta.s) ~ mydata$physa.mean.metas)
summary(DO)
plot(mydata$do.delta.s ~ mydata$physa.mean.metas)
shapiro.test(DO$residuals)



### nutrients in the mesocosms
TN<- read.csv("nutrients.csv", header = TRUE)

#nutrients for week 6 and later (once parasitism could have impact)
nut.late <- TN[25:97,]

### DOC mixed model
library(lattice)
library(lme4)
library(blmeco)
library(car)
library(effects)
M2 <- lmer(DOC ~ physa.mean.metas + (1 | tank), data = nut.late)
summary(M2)
# for reasonable sample sizes, can use Wald inference:
Betas  <- fixef(M2)                  #Get the betas
SE     <-  sqrt(diag(vcov(M2)))      #Get the SEs
pval   <- 2*pnorm(-abs(Betas  / SE)) #Z distribution
Output <- cbind(Betas,SE, pval)
print(Output, digits = 3) # P values for fixed effect
plot(M2)
qqmath(M2)
Anova(M2, type="III")

Cwolf <- lm(mydata$DOC.t ~ mydata$total.mass)
summary(Cwolf)
shapiro.test(Cwolf$residuals)#significance
N <- lm(mydata$TN.t ~ mydata$total.mass)
summary(N)
shapiro.test(N$residuals)
P <- lm(mydata$TDP.t ~ mydata$total.mass)
summary(P)
shapiro.test(P$residuals)

Cwolf <- lm(mydata$DOC.t ~ mydata$dm.per.cm2)
summary(Cwolf)
shapiro.test(Cwolf$residuals)#significance
N <- lm(mydata$TN.t ~ mydata$dm.per.cm2)
summary(N)
shapiro.test(N$residuals)
P <- lm(log(mydata$TDP.t) ~ mydata$dm.per.cm2)
summary(P)
shapiro.test(P$residuals)


#total dissolved nitrogen dynamics
N.mod <- lm(log(TN$TN) ~ TN$physa.mean.metas * TN$week)
summary(N.mod)
plot(N.mod)
shapiro.test(N.mod$residuals)#signifiance

M2 <- lmer(TN ~ physa.mean.metas * week + (1 | tank), data = TN)
summary(M2)
# for reasonable sample sizes, can use Wald inference:
Betas  <- fixef(M2)                  #Get the betas
SE     <-  sqrt(diag(vcov(M2)))      #Get the SEs
pval   <- 2*pnorm(-abs(Betas  / SE)) #Z distribution
Output <- cbind(Betas,SE, pval)
print(Output, digits = 3) # P values for fixed effect
plot(M2)
qqmath(M2)
Anova(M2, type="III")



### TDP model
P.mod <- lm(log(TN$TP) ~ TN$physa.mean.metas * TN$week)
summary(P.mod)
plot(P.mod)
shapiro.test(P.mod$residuals)


### look at stoichmetry
CP.mod <- lm(log(TN$CP) ~ TN$physa.mean.metas * TN$week)
summary(CP.mod)
plot(CP.mod)
shapiro.test(CP.mod$residuals)

CN.mod <- lm(TN$CN ~ TN$physa.mean.metas * TN$week)
summary(CN.mod)
plot(CN.mod)
shapiro.test(CN.mod$residuals)

library(mblm)
model.k = mblm(CN ~ physa.mean.metas,
               data=TN)
summary(model.k)

NP.mod <- lm(log(TN$NP) ~ TN$physa.mean.metas * TN$week)
summary(NP.mod)
plot(NP.mod)
shapiro.test(NP.mod$residuals)

#N:P stoich residuals are not normal. Attempt a box-cox transformation
my.bc <- boxcox(NP.mod, lambda=seq(-0.5, 1, by=0.01))
max(my.bc$y)
my.bc$x[which(my.bc$y==(max(my.bc$y)))]

#suggests raising NP data to -0.04 will be closest to normal.
NP.mod <- lm((TN$NP^0.94) ~ TN$physa.mean.metas * TN$week)
summary(NP.mod)
plot(NP.mod)
shapiro.test(NP.mod$residuals)

#cannot get NP to normal so use nonparametric ancova from sm pacakge; but have to treat treatment as a factor
library(sm)
np.NP.mod <- sm.ancova( TN$week, TN$NP, TN$factor, model='equal')

########################
### Analayses using prevalence as a predictor
########################
physa <- lm(mydata$physa7 ~ mydata$prev)
summary(physa)
shapiro.test(physa$residuals)


### look at snail abundance across the treatment sets
physa <- lm(log(mydata$physa.abun) ~ mydata$prev)
summary(physa)
shapiro.test(physa$residuals)

promen <- lm(mydata$pro.abun ~ mydata$prev)
summary(promen)
shapiro.test(promen$residuals)

total.ab <- lm(log(mydata$total.abun) ~ mydata$prev)
summary(total.ab)
shapiro.test(total.ab$residuals)

phys.pro <- lm(log(mydata$phys.pro) ~ mydata$prev)
summary(phys.pro)
shapiro.test(phys.pro$residuals)

phys.all <- lm((mydata$physa.abun/mydata$total.abun) ~ mydata$prev)
summary(phys.all)
shapiro.test(phys.all$residuals)


### periphyton characteristics
#create quadratic value for infected-snail days
infected2 <- mydata$prev * mydata$prev

#periphyton percent ash-free dry mass
AFDMq <-lm(mydata$X.AFDM ~ mydata$prev)# + infected2)
summary(AFDMq)
plot(AFDMq)
shapiro.test(AFDMq$residuals) #significance

### chlorophyll a
chloro <- lm(mydata$chla.per.cm2 ~ mydata$prev)
summary(chloro)
plot(chloro)
shapiro.test(chloro$residuals)

#periphyton dry mass per square centimeter
periq <-lm(log(mydata$dm.per.cm2) ~ mydata$prev)# + infected2)
summary(periq)
shapiro.test(periq$residuals)
plot(periq)



###Wolfiella biomass
wolf.mod <- lm(mydata$total.wolf ~ mydata$prev)
summary(wolf.mod)
plot(wolf.mod)
shapiro.test(wolf.mod$residuals)

surfveg.mod <- lm(mydata$total.mass ~ mydata$prev)
summary(surfveg.mod)
shapiro.test(surfveg.mod$residuals)

lemna.mod <- lm(log(mydata$total.lemna) ~ mydata$prev)
summary(lemna.mod)
shapiro.test(lemna.mod$residuals)

wolf.lem <- lm(log(mydata$wolf.lem) ~ mydata$prev)
summary(wolf.lem)
shapiro.test(wolf.lem$residuals)

###in situ diel primary production
#week 6
DO <- lm(log(mydata$do.delta.s) ~ mydata$prev)
summary(DO)
plot(mydata$do.delta.s ~ mydata$prev)
shapiro.test(DO$residuals)



### nutrients in the mesocosms
TN<- read.csv("nutrients.csv", header = TRUE)

#nutrients for week 6 and later (once parasitism could have impact)
nut.late <- TN[25:97,]

### DOC mixed model
library(lattice)
library(lme4)
library(blmeco)
library(car)
library(effects)
M2 <- lmer(DOC ~ prev + (1 | tank), data = nut.late)
summary(M2)
# for reasonable sample sizes, can use Wald inference:
Betas  <- fixef(M2)                  #Get the betas
SE     <-  sqrt(diag(vcov(M2)))      #Get the SEs
pval   <- 2*pnorm(-abs(Betas  / SE)) #Z distribution
Output <- cbind(Betas,SE, pval)
print(Output, digits = 3) # P values for fixed effect
plot(M2)
qqmath(M2)
Anova(M2, type="III")

Cwolf <- lm(mydata$DOC.t ~ mydata$total.wolf)
summary(Cwolf)
shapiro.test(Cwolf$residuals)#significance

#total dissolved nitrogen dynamics
N.mod <- lm(log(TN$TN) ~ TN$prev * TN$week)
summary(N.mod)
plot(N.mod)
shapiro.test(N.mod$residuals)




### TDP model
P.mod <- lm(log(TN$TP) ~ TN$prev * TN$week)
summary(P.mod)
plot(P.mod)
shapiro.test(P.mod$residuals)


### look at stoichmetry
CP.mod <- lm(log(TN$CP) ~ TN$prev * TN$week)
summary(CP.mod)
plot(CP.mod)
shapiro.test(CP.mod$residuals)

CN.mod <- lm(TN$CN ~ TN$prev * TN$week)
summary(CN.mod)
plot(CN.mod)
shapiro.test(CN.mod$residuals)


NP.mod <- lm(log(TN$NP) ~ TN$prev * TN$week)
summary(NP.mod)
plot(NP.mod)
shapiro.test(NP.mod$residuals)

#N:P stoich residuals are not normal. Attempt a box-cox transformation
my.bc <- boxcox(NP.mod, lambda=seq(-0.5, 1, by=0.01))
max(my.bc$y)
my.bc$x[which(my.bc$y==(max(my.bc$y)))]

#suggests raising NP data to -0.04 will be closest to normal.
NP.mod <- lm((TN$NP^0.93) ~ TN$prev * TN$week)
summary(NP.mod)
plot(NP.mod)
shapiro.test(NP.mod$residuals)

#cannot get NP to normal so use nonparametric ancova from sm pacakge; but have to treat treatment as a factor
library(sm)
np.NP.mod <- sm.ancova( TN$week, TN$NP, TN$factor, model='equal')




#### Additional analyses within supplementary material ####

### invertebrate community NMDS
library(vegan)
data <- read.csv("community.csv")

#rare species occuring in only one tank or ubiquitous species are removed as these provide little information
#or spurious correlations. Cladocerans are analysed at the genus level.
zoos <- metaMDS(data[,c(6,9:12,14:15,17:23)], distance = "jaccard" , k=2, trymax=100)		## Ordinate in 2 dim, 20 tries for best solution
plot(zoos, type = 't', xlim = c(-1.75,1))
plot(zoos$points)
points(zoos, pch = 16, cex = 2, col=as.numeric(data$col + 2))
stressplot(zoos)	
#prep data for adonis
treat <- as.factor(data$col)
species <- data[,c(6,9:12,14:15,17:23)]
# Do adonis to check for differences
adon <- adonis(species ~ treat, method="jaccard", data=data, control=permControl(strata=treat), permuations=9999)
adon

### TDP model
P.mod <- lm(log(TN$TP) ~ TN$infected.days * TN$week)
summary(P.mod)
plot(P.mod)
shapiro.test(P.mod$residuals)


### look at stoichmetry
CP.mod <- lm(log(TN$CP) ~ TN$infected.days * TN$week)
summary(CP.mod)
plot(CP.mod)
shapiro.test(CP.mod$residuals)

CN.mod <- lm(TN$CN ~ TN$infected.days * TN$week)
summary(CN.mod)
plot(CN.mod)
shapiro.test(CN.mod$residuals)


NP.mod <- lm(TN$NP ~ TN$infected.days * TN$week)
summary(NP.mod)
plot(NP.mod)
shapiro.test(NP.mod$residuals)

#N:P stoich residuals are not normal. Attempt a box-cox transformation
my.bc <- boxcox(NP.mod, lambda=seq(-0.5, 0.5, by=0.01))
max(my.bc$y)
my.bc$x[which(my.bc$y==(max(my.bc$y)))]

#suggests raising NP data to -0.04 will be closest to normal.
NP.mod <- lm((TN$NP^ -0.04) ~ TN$infected.days * TN$week)
summary(NP.mod)
plot(NP.mod)
shapiro.test(NP.mod$residuals)

#cannot get NP to normal so use nonparametric ancova from sm pacakge; but have to treat treatment as a factor
library(sm)
np.NP.mod <- sm.ancova( TN$week, TN$NP, TN$factor, model='equal')


#SEM to parse out influence of density and infection metrics
library(lavaan)
mydata <- read.csv("mesocosm.master.csv", header = TRUE)
#physa.mean.metas and prev are the predictors
model1<- "
	#Regressions
	#dm.per.cm2 ~ 0.27
	dm.per.cm2 ~ total.abun + physa.mean.metas + treatment
	physa.mean.metas ~ treatment
	total.abun ~ physa.mean.metas
	"

sem.fit1<-sem(model1, data=mydata)

summary(sem.fit1, standardized=T)
inspect(sem.fit1, 'r2')

	#dm.per.cm2 ~ total.abun
	#(total.abun ~ physa.mean.metas) * (dm.per.cm2 ~ physa.mean.metas)
	#(total.abun ~ physa.mean.metas) * (dm.per.cm2 ~ physa.mean.metas) + dm.per.cm2 ~ physa.mean.metas


ab <- glm(dm.per.cm2 ~ total.abun, data=mydata, family = Gamma(link = 'inverse'))
summary(ab)
shapiro.test(ab$residuals)
plot(ab)

model2<- "
	#Regressions
	#X.AFDM ~ 28
	X.AFDM ~ prev + treatment + total.abun
	prev ~ treatment
	total.abun ~ prev
	"

sem.fit2<-sem(model2, data=mydata)

summary(sem.fit2, standardized=T)
inspect(sem.fit2, 'r2')
