#Question (4) How do oyster spat densities compare across projects and cultch densities 
#in Apalachicola Bay?

library(readxl)
library(tidyverse)
library(dplyr)
library(Hmisc)
library(MASS)
library(sjPlot)
library(reshape)
library(ggplot2)
library(lubridate)
library(AICcmodavg)
library(ggeffects)
library(cowplot)
library(ggplot2); theme_set(theme_bw(base_size=16) + theme(panel.spacing = grid::unit(0, "lines")))
library(glmmTMB)
library(bbmle)
library(ggeffects)
library(AICcmodavg)
library(emmeans)
library(DHARMa)
library(readr)
library(jtools)

d3 <- read_csv("d3_rock_spat.csv")

d3$SP <- with(d3, interaction(Site, Project, sep = "_", drop = TRUE))

dApalach<-subset(d3,d3$Bay =="Apalachicola")

plot(dApalach$Spat_sum~dApalach$Roundwt)

tab <- with(dApalach,table(SP,Project))
library(Matrix)
ifun <- function(M) {
  M <- as(M, "Matrix")
  image(M, aspect = "fill",
        scales = list(y = list(at = seq(nrow(M)), labels = rownames(M)),
                      x = list(at = seq(ncol(M)), labels = colnames(M), rot = 90)))
}
ifun(tab)



## sum-to-zero contrasts so main effect of Period = unweighted average across Bays
options(contrasts = c("contr.sum", "contr.poly"))

######

plot(dApalach$Spat_sum~dApalach$Roundwt)
m1<-lm(Spat_sum~Roundwt, data=dApalach)
summary(m1)


ggplot(dApalach, aes(x=predict(m1), y= Spat_sum)) +
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  labs(x='Round weight', y='Spat sum', title='')



#########check for overparameterization


#Intercept
tmb0.AB <- glmmTMB(Roundwt ~ (1|SP) + offset(log(Num_quads)),
                   data = dApalach, family="nbinom2") #converge
summary(tmb0.AB)


ggplot(dApalach, aes(x=predict(tmb0.AB), y= Spat_sum)) +
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  labs(x='Round weight', y='Spat sum', title='')


#spat sum

tmb00.AB <- update(tmb0.AB, . ~ . + Spat_sum)
summary(tmb00.AB)


#Period
tmb1.AB <- update(tmb00.AB, . ~ . + Period)
summary(tmb1.AB)

#Period + Project
tmb2.AB <- update(tmb1.AB, . ~ . + Project)
summary(tmb2.AB)


#Period*Project
tmb3.AB <- update(tmb2.AB, . ~ . + Period:Project)
summary(tmb3.AB)


########################

#Project
tmb4.AB <- update(tmb0.AB, . ~ . + Project)
summary(tmb4.AB)



#Period*bay/site (allow period by site across project)
tmb5.AB <- update(tmb3.AB, . ~ . - (1|SP) + (Period|SP))
diagnose(tmb5.AB)  ##  BMB: this is a *singular fit*: correlation of -1
## means this is probably overfitted
VarCorr(tmb5.AB)

summary(tmb5.AB) #tmb 5 has spat sum in it
#and spat sum not significant

#plot coefficients
plot_summs(tmb5.AB)

#tmb 5.ab.xx is spat sum:project
tmb5.AB.xx <- update(tmb3.AB, . ~ . - (1|SP) - Spat_sum + Project:Spat_sum + (Period|SP)) #spat sum*project

diagnose(tmb5.AB.xx)  
VarCorr(tmb5.AB)


summary(tmb5.AB.xx)

(em5.ABxx <- emtrends(tmb5.AB.xx, ~Project, "Spat_sum"))
test(em5.ABxx)


#all + dispersion
tmb6.AB <- update(tmb5.AB, dispformula = ~Project)
summary(tmb6.AB)

#bp note, tmb6.AB is tmb5.AB + adding a unique dispersion parameter for each project. 
#has singular convergence issue goes away with bfgs 

tmb7.AB <- update(tmb3.AB, . ~ .  + (0+Period|SP), dispformula = ~Project)
summary(tmb7.AB)

diagnose(tmb7.AB)
VarCorr(tmb7.AB)

#model selection information

## self-naming list
cand.set2.AB =
  list(tmb0.AB, tmb00.AB, tmb1.AB, tmb2.AB, tmb3.AB, tmb4.AB, tmb5.AB, tmb5.AB.xx, tmb6.AB, tmb7.AB)
modnames2.AB = c("intercept","spat sum", "period", "period + project", "period*project", "project", "all",
                 "all +project:spat_sum","all + dispersion", "project/sp uncorr + disp")
names(cand.set2.AB) <- modnames2.AB


#AIC
aictab(cand.set2.AB, modnames2.AB, second.ord = FALSE) #model selection table with AIC
#AICc
aictab(cand.set2.AB, modnames2.AB, second.ord = TRUE) #model selection table with AICc

#for example, grab tmb5 and fit it w/o spat sum and see if different

tmb5.AB.nospat <- update(tmb3.AB, . ~ . - (1|SP) - Spat_sum + (Period|SP)) #w/o spat sum
AIC(tmb5.AB.nospat,tmb5.AB) #delta AIC about 1.8. so not separable (w/o is slightly lower)


summary(tmb5.AB.nospat)

#it isn't

#5 simpler models in an AIC table
cand.set3.AB =
  list(tmb0.AB, tmb00.AB, tmb1.AB, tmb3.AB, tmb4.AB, tmb5.AB, tmb5.AB.xx, tmb6.AB, tmb5.AB.nospat)
modnames3.AB = c("intercept","spat sum", "period", "period * project",  "project", "all", 
                 "all*project:spat", "all + unqiue disp","all no spat")
names(cand.set3.AB) <- modnames3.AB


#AIC
aictab(cand.set3.AB, modnames3.AB, second.ord = FALSE) #model selection table with AIC
#AICc
aictab(cand.set3.AB, modnames3.AB, second.ord = TRUE) #model selection table with AICc


# ## quantify and test trends by project

#tmb5.AB.xx is without spatsum
#tmb5.AB is with spat sum

(em5.AB <- emtrends(tmb5.AB, ~Project, "Spat_sum"))
test(em5.AB)

summary(tmb5.AB)

#is period sig with spat by project
(em5x.AB <- emtrends(tmb5.AB.xx, ~Project, "Period"))
test(em5x.AB)


