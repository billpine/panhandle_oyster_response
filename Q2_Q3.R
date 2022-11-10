#Questions 2 and 3 in this code. First Q2 then down below Q3


#Question 2 In a focal site (Apalachicola Bay), how do trends in oyster 
#spat (the life stage hypothesized to respond first to restoration) vary 
#among separate restoration projects?

library(ggplot2); theme_set(theme_bw(base_size=16) + theme(panel.spacing = grid::unit(0, "lines")))
library(glmmTMB)
library(bbmle)
library(ggeffects)
library(AICcmodavg)
library(emmeans)
library(DHARMa)
library(readr)
library(jtools)

############################


dp4 <- read.csv("dp4x.csv")

dp4$Site<-as.factor(dp4$Site)
dp4$Project<-as.factor(dp4$Project)


##Apalachicola

head(dp4)

names(dp4)

dp4$SP <- with(dp4, interaction(Site, Project, sep = "_", drop = TRUE))

#check this with a graphic

tab <- with(dp4,table(SP,Project))
library(Matrix)
ifun <- function(M) {
  M <- as(M, "Matrix")
  image(M, aspect = "fill",
        scales = list(y = list(at = seq(nrow(M)), labels = rownames(M)),
                      x = list(at = seq(ncol(M)), labels = colnames(M), rot = 90)))
}
ifun(tab)

###########
##graph

ggplot(dp4, aes(Period, Sum_spat)) + facet_wrap(~Project) + geom_point() + geom_line(aes(group = SP), alpha = 0.5) +
  scale_y_continuous(trans = scales::log1p_trans())

## sum-to-zero contrasts so main effect of Period = unweighted average across Bays
options(contrasts = c("contr.sum", "contr.poly"))


#Intercept
tmb0.AB <- glmmTMB(Sum_spat ~ (1|SP) + offset(log(Num_quads)),
                   data = dp4, family="nbinom2", control=glmmTMBControl(optimizer=optim,
                                                                        optArgs=list(method="BFGS"))) #converge
summary(tmb0.AB)

#Period
tmb1.AB <- update(tmb0.AB, . ~ . + Period, control=glmmTMBControl(optimizer=optim,
                                                                  optArgs=list(method="BFGS")))
summary(tmb1.AB)

#Period + Project
tmb2.AB <- update(tmb1.AB, . ~ . + Project, control=glmmTMBControl(optimizer=optim,
                                                                   optArgs=list(method="BFGS")))
summary(tmb2.AB)

#compare optimizers

tmb2.AB1 <- update(tmb1.AB, . ~ . + Project)

#fixed effects
all.equal(fixef(tmb2.AB1), fixef(tmb2.AB))

#Period*Project
tmb3.AB <- update(tmb2.AB, . ~ . + Period:Project, control=glmmTMBControl(optimizer=optim,
                                                                          optArgs=list(method="BFGS")))
summary(tmb3.AB)

###just to compare optimizers
#default nlminb
tmb3.AB1 <- update(tmb2.AB, . ~ . + Period:Project)
#bfgs
tmb3.AB2 <- update(tmb2.AB, . ~ . + Period:Project, control=glmmTMBControl(optimizer=optim,
                                                                           optArgs=list(method="BFGS")))
#fixed effects
all.equal(fixef(tmb3.AB1), fixef(tmb3.AB2))
#yes all equal


#plot coefficients
plot_summs(tmb3.AB2)

########################

#Project
tmb4.AB <- update(tmb0.AB, . ~ . + Project, control=glmmTMBControl(optimizer=optim,
                                                                   optArgs=list(method="BFGS")))
summary(tmb4.AB)
#Period*bay/site (allow period by site across project)
tmb5.AB <- update(tmb3.AB, . ~ . - (1|SP) + (Period|SP),control=glmmTMBControl(optimizer=optim,
                                                                               optArgs=list(method="BFGS")))
diagnose(tmb5.AB)  ##  BMB: this is a *singular fit*: correlation of -1
## means this is probably overfitted
VarCorr(tmb5.AB)

summary(tmb5.AB)

#plot coefficients
plot_summs(tmb5.AB)


#this is just tmb5.AB but adds a unique dispersion term for each project 

#all + dispersion
tmb6.AB <- update(tmb5.AB, dispformula = ~Project,control=glmmTMBControl(optimizer=optim,
                                                                         optArgs=list(method="BFGS")))
summary(tmb6.AB)

tmb7.AB <- update(tmb3.AB, . ~ .  + (0+Period|SP), dispformula = ~Project, control=glmmTMBControl(optimizer=optim,
                                                                                                  optArgs=list(method="BFGS")))
summary(tmb7.AB)

#model selection information

## self-naming list
cand.set2.AB =
  list(tmb0.AB,tmb1.AB,tmb2.AB,tmb3.AB,tmb4.AB, tmb5.AB, tmb6.AB, tmb7.AB)
modnames2.AB = c("intercept", "period", "period + project", "period*project", "project", "all",
                 "all + dispersion", "project/sp uncorr + disp")
names(cand.set2.AB) <- modnames2.AB


#AIC
aictab(cand.set2.AB, modnames2.AB, second.ord = FALSE) #model selection table with AIC
#AICc
aictab(cand.set2.AB, modnames2.AB, second.ord = TRUE) #model selection table with AICc

#AIC(c) table of all models
AICctab(cand.set2.AB, weights=TRUE)


#plot coefficients
plot_summs(tmb5.AB, tmb6.AB)


## quantify and test trends by project
(em1.AB <- emtrends(tmb5.AB, ~Project, "Period"))
test(em1.AB)

#compare to the project*period model
(em3.AB <- emtrends(tmb3.AB, ~Project, "Period"))
test(em3.AB)

#predict by project or period

ggpredict(tmb5.AB)
test.nfwf1 = ggpredict(tmb5.AB, terms = c("Period[9]", "Project[NFWF-1]", "Num_quads[1]"), type = c('fe')) 
test.nrda4044 = ggpredict(tmb5.AB, terms = c("Period[13]", "Project[NRDA-4044]", "Num_quads[1]"), type = c('fe')) 
test.nrda5077 = ggpredict(tmb5.AB, terms = c("Period[12]", "Project[GEBF-5007]", "Num_quads[1]"), type = c('fe')) 
test.fwc2021 = ggpredict(tmb5.AB, terms = c("Period[15]", "Project[NFWF-2021]", "Num_quads[1]"), type = c('fe')) 

#########
#check autocorrelation of best
#
#########

#check autocorrelation issues
res1 <- simulateResiduals(tmb5.AB)
plot(res1)
agg.res1 = recalculateResiduals(res1,group=dp4$Period)
time = unique(dp4$Period)
plot(time,agg.res1$scaledResiduals,pch=16)

#######
###
#Question 3 Are oyster spat counts in Apalachicola Bay
#associated with freshwater discharge? 


#now tmb5 w. discharge
tmb5.12k <- update(tmb5.AB, . ~ . + Lowdays_12, control=glmmTMBControl(optimizer=optim,
                                                                       optArgs=list(method="BFGS")))
summary(tmb5.12k)

tmb5.lag12 <- update(tmb5.AB, . ~ . + lag1, control=glmmTMBControl(optimizer=optim,
                                                                   optArgs=list(method="BFGS")))
summary(tmb5.lag12)
tmb5.6k <- update(tmb5.AB, . ~ . + Lowdays_6, control=glmmTMBControl(optimizer=optim,
                                                                     optArgs=list(method="BFGS")))
tmb5.lag6 <- update(tmb5.AB, . ~ . + lag1_6, control=glmmTMBControl(optimizer=optim,
                                                                    optArgs=list(method="BFGS")))

cand.set3.1.AB =
  list(tmb5.AB,tmb5.12k,tmb5.lag12, tmb5.6k, tmb5.lag6)
modnames3.1.AB = c("full", "full_low12k", "full_12k_lag", "full_low6k", "full_6k_lag")
names(cand.set3.1.AB) <- modnames3.1.AB

#AICc
aictab(cand.set3.1.AB, modnames3.1.AB, second.ord = TRUE) #model selection table with AICc

#no improvement with the different discharge metrics

###############

