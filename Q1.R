#Question 1 How do temporal trends in oyster counts vary among the 
#three depressed bays where restoration has taken place 
#(Pensacola, St. Andrew, and Apalachicola bays)?

library(ggplot2); theme_set(theme_bw(base_size=16) + theme(panel.spacing = grid::unit(0, "lines")))
library(glmmTMB)
library(bbmle)
library(ggeffects)
library(AICcmodavg)
library(emmeans)
library(DHARMa)
library(readr)
library(jtools)

d5 <- read_csv("d5.csv")
d5$Bay <- as.factor(d5$Bay)
str(d5)

tab <- with(d5,table(Site,Bay))
library(Matrix)
ifun <- function(M) {
    M <- as(M, "Matrix")
    image(M, aspect = "fill",
          scales = list(y = list(at = seq(nrow(M)), labels = rownames(M)),
                        x = list(at = seq(ncol(M)), labels = colnames(M), rot = 90)))
}
ifun(tab)

options(contrasts = c("contr.sum", "contr.poly"))

#Intercept
tmb0 <- glmmTMB(Sum_spat ~ (1|Site) + offset(log(Num_quads)),
                data = d5, family="nbinom2") #converge
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

#tmb6
#all + dispersion

tmb6 <- update(tmb5, dispformula = ~Bay)
summary(tmb6)

ggplot(d5, aes(Period, Sum_spat)) + facet_wrap(~Bay) + geom_point() + geom_line(aes(group = Site), alpha = 0.5) +
  scale_y_continuous(trans = scales::log1p_trans())

cand.set2 = list(tmb0,tmb1,tmb2,tmb3,tmb4, tmb5, tmb6)
modnames2 = c("intercept", "period", "period + bay", "period*bay", "bay", "period*bay/site",
              "all + dispersion")
names(cand.set2) <- modnames2


#AIC(c) table of all models
AICctab(cand.set2, weights=TRUE)

aictab(cand.set2, modnames2, second.ord = FALSE) #model selection table with AIC
aictab(cand.set2, modnames2, second.ord = TRUE) #model selection table with AICc


#check autocorrelation issues
res1 <- simulateResiduals(tmb5)
plot(res1)
agg.res1 = recalculateResiduals(res1,group=d5$Period)
time = unique(dp4$Period)
plot(time,agg.res1$scaledResiduals,pch=16)

## quantify and test trends by bay

(em1 <- emtrends(tmb5, ~Bay, "Period"))
test(em1)

#marginal means for tmb5
mm.tmb5<-emmeans(tmb5, ~Bay, "Period")

library(ggeffects)
ggpredict(tmb5)
pred_tmb5 <- ggpredict(tmb5, c("Period[15]", "Bay","Num_quads[1]"))



############################
############################

#now for seed and legal (using best model from above)
tmb5.1_seed <- glmmTMB(Sum_seed ~ Period + Bay + (Period|Site) + Period:Bay + offset(log(Num_quads)),
                  data = d5, family="nbinom2") #converge
summary(tmb5.1_seed)

(em1_seed <- emtrends(tmb5.1_seed, ~Bay, "Period"))
test(em1_seed)

pred_tmb5.1_seed <- ggpredict(tmb5.1_seed, c("Period[15]", "Bay","Num_quads[1]"))

##
tmb5.1_legal <- glmmTMB(Sum_legal ~ Period + Bay + (Period|Site) + Period:Bay + offset(log(Num_quads)),
                       data = d5, family="nbinom2",control=glmmTMBControl(optimizer=optim,
                                                                          optArgs=list(method="BFGS"))) #converge
summary(tmb5.1_legal)

(em1_legal <- emtrends(tmb5.1_legal, ~Bay, "Period"))
test(em1_legal)

pred_tmb5.1_legal <- ggpredict(tmb5.1_legal, c("Period[15]", "Bay","Num_quads[1]"))


############################
