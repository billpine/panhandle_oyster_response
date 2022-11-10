#Q5 How well do different types and densities of restoration-sourced 
#cultch persist following deployment in Pensacola, St. Andrew, and 
#Apalachicola bays? 

d3 <- read_csv("d3_rock.csv")

#######
#######


as.integer(d3$Roundwt)


d3$SP <- with(d3, interaction(Site, Project, sep = "_", drop = TRUE))

tab <- with(d3,table(Site,Bay))
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

ggplot(d3, aes(Period, Roundwt)) + facet_wrap(~Bay) + geom_point() + geom_line(aes(group = Site), alpha = 0.5) +
  scale_y_continuous(trans = scales::log1p_trans())

#Intercept
tmb0 <- glmmTMB(Roundwt ~ (1|Site) + offset(log(Num_quads)),
                data = d3, family="nbinom2") #converge
summary(tmb0)

#Period
tmb1 <- update(tmb0, . ~ . + Period)
summary(tmb1)

#Period + Bay
tmb2 <- update(tmb1, . ~ . + Bay)
summary(tmb2)

#Period*Bay
tmb3 <- update(tmb2, . ~ . + Period:Bay)
summary(tmb3)

#Bay
tmb4 <- update(tmb0, . ~ . + Bay)
summary(tmb4)


#Period*bay/site (allow period by site across bay)
tmb5 <- update(tmb3, . ~ . - (1|Site) + (Period|Site))

summary(tmb5)

# fails with bfgs as well
tmb5.x <- update(tmb3, . ~ . - (1|Site) + (Period|Site), control=glmmTMBControl(optimizer=optim,
                                                                                optArgs=list(method="BFGS")))


tmb6 <- update(tmb5, dispformula = ~Bay)
                                                              optArgs=list(method="BFGS")))
summary(tmb6)

#drop tmb5 & 6 bc of convergence (w/ either optimizer)
cand.set2 = list(tmb0,tmb1,tmb2,tmb3,tmb4)
modnames2 = c("intercept", "period", "period + bay", "period*bay", "bay")
names(cand.set2) <- modnames2

#model selection information

aictab(cand.set2, modnames2, second.ord = FALSE) #model selection table with AIC
aictab(cand.set2, modnames2, second.ord = TRUE) #model selection table with AICc

(em1 <- emtrends(tmb5, ~Bay, "Period"))
test(em1)

ggpredict(tmb5)

pred_tmb5<- ggpredict(tmb5, c("Period[15]", "Bay","Num_quads[1]"))

#########
#########
#Just working with Apalachicola data
#########
#########

dApalach<-subset(d3,d3$Bay =="Apalachicola")

dApalach$SP <- with(dApalach, interaction(Site, Project, sep = "_", drop = TRUE))

tab <- with(dApalach,table(SP,Project))
library(Matrix)
ifun <- function(M) {
  M <- as(M, "Matrix")
  image(M, aspect = "fill",
        scales = list(y = list(at = seq(nrow(M)), labels = rownames(M)),
                      x = list(at = seq(ncol(M)), labels = colnames(M), rot = 90)))
}
ifun(tab)

names(dApalach)

ggplot(dApalach, aes(Period, Roundwt)) + facet_wrap(~Project) + geom_point() + geom_line(aes(group = SP), alpha = 0.5) +
  scale_y_continuous(trans = scales::log1p_trans())

#Intercept
tmb0.AB <- glmmTMB(Roundwt ~ (1|SP) + offset(log(Num_quads)),
                   data = dApalach, family="nbinom2") #converge
summary(tmb0.AB)


#Period
tmb1.AB <- update(tmb0.AB, . ~ . + Period)
summary(tmb1.AB)

#Period + Project
tmb2.AB <- update(tmb1.AB, . ~ . + Project)
summary(tmb2.AB)

#Period*Project
tmb3.AB <- update(tmb2.AB, . ~ . + Period:Project)
summary(tmb3.AB)

#Project
tmb4.AB <- update(tmb0.AB, . ~ . + Project)
summary(tmb4.AB)

#Period*Project/site (allow period by site across Project)
tmb5.AB <- update(tmb3.AB, . ~ . - (1|SP) + (Period|SP))
diagnose(tmb5.AB)  
VarCorr(tmb5.AB)

summary(tmb5.AB)

#plot coefficients
plot_summs(tmb5.AB)

#all + dispersion
tmb6.AB <- update(tmb5.AB, dispformula = ~Project)
summary(tmb6.AB)

tmb7.AB <- update(tmb3.AB, . ~ .  + (0+Period|SP), dispformula = ~Project)
summary(tmb7.AB)

diagnose(tmb7.AB)
VarCorr(tmb7.AB)

#model selection information

## self-naming list
cand.set2 =
  list(tmb0.AB, tmb1.AB, tmb2.AB, tmb3.AB, tmb4.AB, tmb5.AB, tmb6.AB, tmb7.AB)
modnames2 = c("intercept", "period", "period + Project", "period*Project", "Project", "all",
              "all + dispersion", "Project/sp uncorr + disp")
names(cand.set2) <- modnames2


#AIC
aictab(cand.set2, modnames2, second.ord = FALSE) #model selection table with AIC
#AICc
aictab(cand.set2, modnames2, second.ord = TRUE) #model selection table with AICc


(em1 <- emtrends(tmb3.AB, ~Project, "Period"))
test(em1)

#plot coefficients
plot_summs(tmb3.AB)


test.nfwf1 = ggpredict(tmb3.AB, terms = c("Period[9]", "Project[NFWF-1]", "Num_quads[1]"), type = c('fe')) 
test.nrda4044 = ggpredict(tmb3.AB, terms = c("Period[13]", "Project[NRDA-4044]", "Num_quads[1]"), type = c('fe')) 
test.nrda5077 = ggpredict(tmb3.AB, terms = c("Period[12]", "Project[GEBF-5007]", "Num_quads[1]"), type = c('fe')) 
test.fwc2021 = ggpredict(tmb3.AB, terms = c("Period[15]", "Project[FWC-2021]", "Num_quads[1]"), type = c('fe')) 
