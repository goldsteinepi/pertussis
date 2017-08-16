#######################
# Analysis vaccination coverage
# Citation: Goldstein ND, Newbern EC, Evans AA, Drezner K, Welles SL. Choice of measure of vaccination and estimates of risk of pediatric pertussis.
# 7/22/14 -- Neal Goldstein
#######################


### FUNCTIONS ###

library(psych) #describe, describeBy
library(gmodels) #CrossTable, contrast
library(car) #vif
#library(mice) #multiple imputation
library(boot) #bootstrapping

#boot strap 95% CI for AIC, BIC, and adjusted McFadden's pseudo r sq
#called from boot with arguments dataframe, index of samples to pull, and model to run
#note on pseudo r sq, may not be good to present for logistic, see: http://stats.stackexchange.com/questions/3559/which-pseudo-r2-measure-is-the-one-to-report-for-logistic-regression-cox-s
bootAIC = function(data, index, model)
{
  newdata = data[index,]
  if (model==1) {
    model = glm(case~as.factor(pertussis_vax_utd)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=newdata,family=binomial(link="logit"))
    aic = AIC(model)
    bic = BIC(model)
    r2 = 1 - ((model$deviance-3)/model$null.deviance)
    return(c(aic,bic,r2))
  } else if (model==2) {
    model = glm(case~as.factor(pertussis_vax_utd_delayed)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=newdata,family=binomial(link="logit"))
    aic = AIC(model)
    bic = BIC(model)
    r2 = 1 - ((model$deviance-3)/model$null.deviance)
    return(c(aic,bic,r2))
  } else if (model==3) {
    model = glm(case~as.factor(pertussis_vax)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=newdata,family=binomial(link="logit"))
    aic = AIC(model)
    bic = BIC(model)
    r2 = 1 - ((model$deviance-3)/model$null.deviance)
    return(c(aic,bic,r2))
  } else if (model==4) {
    model = glm(case~pertussis_vax+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=newdata,family=binomial(link="logit"))
    aic = AIC(model)
    bic = BIC(model)
    r2 = 1 - ((model$deviance-3)/model$null.deviance)
    return(c(aic,bic,r2))
  } else if (model==5) {
    model = glm(case~as.factor(pertussis_vax_utd)+as.factor(pertussis_vax)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=newdata,family=binomial(link="logit"))
    aic = AIC(model)
    bic = BIC(model)
    r2 = 1 - ((model$deviance-4)/model$null.deviance)
    return(c(aic,bic,r2))
  } else if (model==6) {
    model = glm(case~as.factor(pertussis_vax_utd_delayed)+as.factor(pertussis_vax)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=newdata,family=binomial(link="logit"))
    aic = AIC(model)
    bic = BIC(model)
    r2 = 1 - ((model$deviance-4)/model$null.deviance)
    return(c(aic,bic,r2))
  } else {
    stop("Unknown model")
  }
}


### READ DATA ###

load("study_population_analysis.RData")


### SUBSET DATA ###

analysis = NA
matchid = studypop_final$matchid[1]
case_asis = studypop_final$case_asis[1]
case_age = studypop_final$age_weeks[1]

#only include cases that are classified as confirmed or probable per the AS IS definition
#and can be classified as UTD/NUTD (min age 13 weeks)
for (i in 1:nrow(studypop_final))
{
  cat("\n\n**************** Observation: ",i," ****************\n",sep="")
  
  if ((case_asis != 0) && (case_age>=13))
    analysis = rbind(analysis, studypop_final[i, ])
  
  if ((matchid != studypop_final$matchid[i+1]) && (!is.na(studypop_final$matchid[i+1])))
  {
    matchid = studypop_final$matchid[i+1]
    case_asis = studypop_final$case_asis[i+1]
    case_age = studypop_final$age_weeks[i+1]
  }
}
rm(i,case_asis,matchid,case_age)
analysis = analysis[-1,]
rownames(analysis) = NULL


### RECODE ###

analysis$mother_parous = ifelse(analysis$mother_parity>=1, 1, 0)
analysis$pertussis_vax = ifelse(analysis$pertussis_vax>=5, 5, analysis$pertussis_vax)
analysis$age_1yr = ifelse(analysis$age==0,0,1)

#create a time since last dose measure
analysis$pertussis_vax_last = NA
for (i in 1:nrow(analysis))
{
  if (analysis$pertussis_vax[i]==0) {
    analysis$pertussis_vax_last[i] = NA
  } else if (analysis$pertussis_vax[i]==1) {
    analysis$pertussis_vax_last[i] = analysis$age_weeks[i] - analysis$pertussis_vax1_age[i]
  } else if (analysis$pertussis_vax[i]==2) {
    analysis$pertussis_vax_last[i] = analysis$age_weeks[i] - analysis$pertussis_vax2_age[i]
  } else if (analysis$pertussis_vax[i]==3) {
    analysis$pertussis_vax_last[i] = analysis$age_weeks[i] - analysis$pertussis_vax3_age[i]
  } else if (analysis$pertussis_vax[i]==4) {
    analysis$pertussis_vax_last[i] = analysis$age_weeks[i] - analysis$pertussis_vax4_age[i]
  } else if (analysis$pertussis_vax[i]==5) {
    analysis$pertussis_vax_last[i] = analysis$age_weeks[i] - analysis$pertussis_vax5_age[i]
  } else if (analysis$pertussis_vax[i]==6) {
    analysis$pertussis_vax_last[i] = analysis$age_weeks[i] - analysis$pertussis_vax6_age[i]
  }
}
rm(i)

#log transform time since last dose for normality
analysis$pertussis_vax_last_log = log(analysis$pertussis_vax_last)


### ANALYSIS: Study population ###

describe(analysis$total_vax)
describe(analysis$pertussis_vax)
CrossTable(analysis$pertussis_vax==0)
CrossTable(analysis$pertussis_vax==1)
CrossTable(analysis$pertussis_vax==2)
CrossTable(analysis$pertussis_vax==3)
CrossTable(analysis$pertussis_vax==4)
CrossTable(analysis$pertussis_vax==5)
CrossTable(analysis$pertussis_vax_utd)
CrossTable(analysis$pertussis_vax_utd_delayed)
describe(analysis$pertussis_vax_last)

CrossTable(analysis$age_1yr)
CrossTable(analysis$gender)
CrossTable(analysis$race)
CrossTable(analysis$ethnicity)
describe(analysis$mother_age)
CrossTable(analysis$mother_married)
CrossTable(analysis$mother_foreign_born)
CrossTable(analysis$mother_education)
CrossTable(analysis$mother_insurance)
CrossTable(analysis$mother_parous)


### ANALYSIS: Comparison of cases and controls ###

describeBy(analysis$total_vax, analysis$case); wilcox.test(analysis$total_vax ~ analysis$case)
describeBy(analysis$pertussis_vax, analysis$case); wilcox.test(analysis$pertussis_vax ~ analysis$case)
CrossTable(analysis$pertussis_vax==0, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$pertussis_vax==1, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$pertussis_vax==2, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$pertussis_vax==3, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$pertussis_vax==4, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$pertussis_vax==5, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$pertussis_vax_utd, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$pertussis_vax_utd_delayed, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(analysis$pertussis_vax_last, analysis$case); wilcox.test(analysis$pertussis_vax_last ~ analysis$case)

CrossTable(analysis$age_1yr, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$gender, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$race, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$ethnicity, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(analysis$mother_age, analysis$case); t.test(analysis$mother_age ~ analysis$case)
CrossTable(analysis$mother_married, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$mother_foreign_born, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$mother_education, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$mother_insurance, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$mother_parous, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)


### ANALYSIS: Comparison of UTD ###

CrossTable(analysis$age_1yr, analysis$pertussis_vax_utd, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$gender, analysis$pertussis_vax_utd, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$race, analysis$pertussis_vax_utd, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$ethnicity, analysis$pertussis_vax_utd, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(analysis$mother_age, analysis$pertussis_vax_utd); t.test(analysis$mother_age ~ analysis$pertussis_vax_utd)
CrossTable(analysis$mother_married, analysis$pertussis_vax_utd, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$mother_foreign_born, analysis$pertussis_vax_utd, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$mother_education, analysis$pertussis_vax_utd, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$mother_insurance, analysis$pertussis_vax_utd, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$mother_parous, analysis$pertussis_vax_utd, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)


### ANALYSIS: Comparison by missingness ###

#missing marker
analysis$missing = ifelse(is.na(analysis$age_1yr) | is.na(analysis$race) | is.na(analysis$mother_parous) | is.na(analysis$pertussis_vax_utd), 1, 0)

CrossTable(analysis$age_1yr, analysis$missing, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$gender, analysis$missing, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$race, analysis$missing, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$ethnicity, analysis$missing, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(analysis$mother_age, analysis$missing); t.test(analysis$mother_age ~ analysis$missing)
CrossTable(analysis$mother_married, analysis$missing, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$mother_foreign_born, analysis$missing, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$mother_education, analysis$missing, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$mother_insurance, analysis$missing, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$mother_parous, analysis$missing, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)


### ANALYSIS: Potential confounder selection ###

#age is always included as it was a frequency matched variable

#crude (referent is UTD as that is the most widely used metric)
summary(glm(case~as.factor(pertussis_vax_utd),data=analysis,family=binomial(link="logit")))

#check for associations with case status, p<0.20
summary(glm(case~as.factor(gender),data=analysis,family=binomial(link="logit")))
summary(glm(case~as.factor(race),data=analysis,family=binomial(link="logit")))
summary(glm(case~as.factor(ethnicity),data=analysis,family=binomial(link="logit")))
#summary(glm(case~mother_age,data=analysis,family=binomial(link="logit")))
#summary(glm(case~as.factor(mother_married),data=analysis,family=binomial(link="logit")))
summary(glm(case~as.factor(mother_foreign_born),data=analysis,family=binomial(link="logit")))
#summary(glm(case~as.factor(mother_education),data=analysis,family=binomial(link="logit")))
#summary(glm(case~as.factor(mother_insurance),data=analysis,family=binomial(link="logit")))
summary(glm(case~as.factor(mother_parous),data=analysis,family=binomial(link="logit")))

#check for associations with utd status (contingent on case status), p<0.20
#summary(glm(pertussis_vax_utd~as.factor(gender),data=analysis,family=binomial(link="logit")))
summary(glm(pertussis_vax_utd~as.factor(race),data=analysis,family=binomial(link="logit")))
#summary(glm(pertussis_vax_utd~as.factor(ethnicity),data=analysis,family=binomial(link="logit")))
#summary(glm(pertussis_vax_utd~mother_age,data=analysis,family=binomial(link="logit")))
#summary(glm(pertussis_vax_utd~as.factor(mother_married),data=analysis,family=binomial(link="logit")))
#summary(glm(pertussis_vax_utd~as.factor(mother_foreign_born),data=analysis,family=binomial(link="logit")))
#summary(glm(pertussis_vax_utd~as.factor(mother_education),data=analysis,family=binomial(link="logit")))
#summary(glm(pertussis_vax_utd~as.factor(mother_insurance),data=analysis,family=binomial(link="logit")))
summary(glm(pertussis_vax_utd~as.factor(mother_parous),data=analysis,family=binomial(link="logit")))

#full model with potential confounders
summary(glm(case~as.factor(pertussis_vax_utd)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))

#check for multicollinearity, VIF>=10
vif(glm(case~as.factor(pertussis_vax_utd)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))

#backward remove nonsignificant vars and check for change in estimate >10%
summary(glm(case~as.factor(pertussis_vax_utd)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))

#test for combining vaccination predictors
model1 = glm(case~as.factor(pertussis_vax_utd)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit"))
model2 = glm(case~as.factor(pertussis_vax_utd)+as.factor(pertussis_vax)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit"))
anova(model1,model2,test="Chisq")
rm(model1,model2)
# 
# model1 = glm(case~as.factor(pertussis_vax_utd_delayed)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit"))
# model2 = glm(case~as.factor(pertussis_vax_utd_delayed)+as.factor(pertussis_vax)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit"))
# anova(model1,model2,test="Chisq")
# rm(model1,model2)

#tests for interaction of age*UTD; no evidence
modelNoInteraction = glm(case~as.factor(pertussis_vax_utd)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit"))
modelInteraction = glm(case~as.factor(pertussis_vax_utd)*as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit"))
anova(modelInteraction,modelNoInteraction,test="Chisq")
rm(modelInteraction,modelNoInteraction)

#tests for interaction of race*UTD; no evidence
modelNoInteraction = glm(case~as.factor(pertussis_vax_utd)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit"))
modelInteraction = glm(case~as.factor(pertussis_vax_utd)*as.factor(race)+as.factor(age_1yr)+as.factor(mother_parous),data=analysis,family=binomial(link="logit"))
anova(modelInteraction,modelNoInteraction,test="Chisq")
rm(modelInteraction,modelNoInteraction)

#tests for interaction of parity*UTD; no evidence
modelNoInteraction = glm(case~as.factor(pertussis_vax_utd)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit"))
modelInteraction = glm(case~as.factor(pertussis_vax_utd)*as.factor(mother_parous)+as.factor(age_1yr)+as.factor(race),data=analysis,family=binomial(link="logit"))
anova(modelInteraction,modelNoInteraction,test="Chisq")
rm(modelInteraction,modelNoInteraction)

#test for interaction of UTD*doses; no evidence
modelNoInteraction = glm(case~as.factor(pertussis_vax_utd)+as.factor(pertussis_vax)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit"))
modelInteraction = glm(case~as.factor(pertussis_vax_utd)*as.factor(pertussis_vax)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit"))
anova(modelInteraction,modelNoInteraction,test="Chisq")
rm(modelInteraction,modelNoInteraction)

#final models
summary(glm(case~as.factor(pertussis_vax_utd)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))
summary(glm(case~as.factor(pertussis_vax_utd_delayed)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))
summary(glm(case~as.factor(pertussis_vax)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))
summary(glm(case~pertussis_vax+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))
summary(glm(case~as.factor(pertussis_vax_utd)+as.factor(pertussis_vax)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))
#summary(glm(case~as.factor(pertussis_vax_utd_delayed)+as.factor(pertussis_vax)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))
#summary(glm(case~pertussis_vax_last_log+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))

#test the contrast to see if delayed UTD is truly different than on time UTD
#format for contrast test is for each possible category, specify vector: 0 to ignore,-1 and 1 to compare those two categories
modelContrast = glm(case~as.factor(pertussis_vax_utd_delayed)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit"))
fit.contrast(modelContrast,"as.factor(pertussis_vax_utd_delayed)",c(0,-1,1))


### ANALYSIS: 4th DOSE ###

#summary by 4th dose in children 2 and younger
summary(glm(case~as.factor(pertussis_vax)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=subset(analysis,age_weeks<=104),family=binomial(link="logit")))

#store vaccine effectiveness estimates and CI
ve=NULL
ve_low=NULL
ve_high=NULL

#loop through age 2 yr (104 weeks) through max age
for (i in 104:max(analysis$age_weeks))
{
  cat("\n\n**************** Observation: ",i," ****************\n",sep="")
  
  #subset analysis
  dose4 = analysis[analysis$age_weeks<=i, ]
  model = glm(case~as.factor(pertussis_vax)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=dose4,family=binomial(link="logit"))
  ve = c(ve, 1 - exp(coef(model)[["as.factor(pertussis_vax)4"]]))
  ve_ci = confint(model)
  ve_low = c(ve_low, 1 - exp(ve_ci["as.factor(pertussis_vax)4",1]))
  ve_high = c(ve_high, 1 - exp(ve_ci["as.factor(pertussis_vax)4",2]))
}
rm(i,dose4,model,ve_ci)

#plot
plot(ve, xlim=c(0,250), ylim=c(0.7,1.0), type="n", axes=F, xlab="Maximum age in study (yrs)", ylab="VE (1-OR)", main="Change in VE for 4th DTaP dose by age of population")
axis(side=1, las=1, at=c(0,52,104,157,209,365), labels=c("2","3","4","5","6","7"))
axis(side=2) 

#smooth out VE lines and put confidence intervals
lines(lowess(ve), lty=1, lwd=1.5)
lines(lowess(ve_low), lty=2)
lines(lowess(ve_high), lty=2)


### COMPLETE CASE ###

#store AIC
results = data.frame(Model=NULL, AIC=NULL, LowerAIC=NULL, UpperAIC=NULL, BIC=NULL, LowerBIC=NULL, UpperBIC=NULL)

#UTD

#regression modeling
summary(glm(case~as.factor(pertussis_vax_utd)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))
confint(glm(case~as.factor(pertussis_vax_utd)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))

#boot strap the 95% CI and mean for AIC & BIC
analysis_boot = boot(analysis, bootAIC, R=1000, model=1)
#mean(analysis_mi_boot$t[,3])
results = rbind(results, data.frame(Model="UTD", AIC=mean(analysis_boot$t[,1]), LowerAIC=boot.ci(analysis_boot, type="norm", index=1)$normal[,2], UpperAIC=boot.ci(analysis_boot, type="norm", index=1)$normal[,3], BIC=mean(analysis_boot$t[,2]), LowerBIC=boot.ci(analysis_boot, type="norm", index=2)$normal[,2], UpperBIC=boot.ci(analysis_boot, type="norm", index=2)$normal[,3]))

#Delayed UTD

#regression modeling
summary(glm(case~as.factor(pertussis_vax_utd_delayed)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))
confint(glm(case~as.factor(pertussis_vax_utd_delayed)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))

#boot strap the 95% CI and mean for AIC & BIC
analysis_boot = boot(analysis, bootAIC, R=1000, model=2)
#mean(analysis_boot$t[,3])
results = rbind(results, data.frame(Model="Delayed UTD", AIC=mean(analysis_boot$t[,1]), LowerAIC=boot.ci(analysis_boot, type="norm", index=1)$normal[,2], UpperAIC=boot.ci(analysis_boot, type="norm", index=1)$normal[,3], BIC=mean(analysis_boot$t[,2]), LowerBIC=boot.ci(analysis_boot, type="norm", index=2)$normal[,2], UpperBIC=boot.ci(analysis_boot, type="norm", index=2)$normal[,3]))

#Categorical shots

#regression modeling
summary(glm(case~as.factor(pertussis_vax)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))
confint(glm(case~as.factor(pertussis_vax)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))

#boot strap the 95% CI and mean for AIC & BIC
analysis_boot = boot(analysis, bootAIC, R=1000, model=3)
#mean(analysis_boot$t[,3])
results = rbind(results, data.frame(Model="Categorical doses", AIC=mean(analysis_boot$t[,1]), LowerAIC=boot.ci(analysis_boot, type="norm", index=1)$normal[,2], UpperAIC=boot.ci(analysis_boot, type="norm", index=1)$normal[,3], BIC=mean(analysis_boot$t[,2]), LowerBIC=boot.ci(analysis_boot, type="norm", index=2)$normal[,2], UpperBIC=boot.ci(analysis_boot, type="norm", index=2)$normal[,3]))

#Continuous shots

#regression modeling
summary(glm(case~pertussis_vax+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))
confint(glm(case~pertussis_vax+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))

#boot strap the 95% CI and mean for AIC & BIC
analysis_boot = boot(analysis, bootAIC, R=1000, model=4)
#mean(analysis_boot$t[,3])
results = rbind(results, data.frame(Model="Continuous doses", AIC=mean(analysis_boot$t[,1]), LowerAIC=boot.ci(analysis_boot, type="norm", index=1)$normal[,2], UpperAIC=boot.ci(analysis_boot, type="norm", index=1)$normal[,3], BIC=mean(analysis_boot$t[,2]), LowerBIC=boot.ci(analysis_boot, type="norm", index=2)$normal[,2], UpperBIC=boot.ci(analysis_boot, type="norm", index=2)$normal[,3]))

#UTD + Categorical shots

#regression modeling
summary(glm(case~as.factor(pertussis_vax_utd)+as.factor(pertussis_vax)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))
confint(glm(case~as.factor(pertussis_vax_utd)+as.factor(pertussis_vax)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))

#boot strap the 95% CI and mean for AIC & BIC
analysis_boot = boot(analysis, bootAIC, R=1000, model=5)
#mean(analysis_boot$t[,3])
results = rbind(results, data.frame(Model="UTD +\nCategorical doses", AIC=mean(analysis_boot$t[,1]), LowerAIC=boot.ci(analysis_boot, type="norm", index=1)$normal[,2], UpperAIC=boot.ci(analysis_boot, type="norm", index=1)$normal[,3], BIC=mean(analysis_boot$t[,2]), LowerBIC=boot.ci(analysis_boot, type="norm", index=2)$normal[,2], UpperBIC=boot.ci(analysis_boot, type="norm", index=2)$normal[,3]))

# #Delayed UTD + Categorical shots
# 
# #regression modeling
# summary(glm(case~as.factor(pertussis_vax_utd_delayed)+as.factor(pertussis_vax)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))
# confint(glm(case~as.factor(pertussis_vax_utd_delayed)+as.factor(pertussis_vax)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=analysis,family=binomial(link="logit")))
# 
# #boot strap the 95% CI and mean for AIC & BIC
# analysis_boot = boot(analysis, bootAIC, R=1000, model=6)
# #mean(analysis_boot$t[,3])
# results = rbind(results, data.frame(Model="Delayed UTD +\nCategorical doses", AIC=mean(analysis_boot$t[,1]), LowerAIC=boot.ci(analysis_boot, type="norm", index=1)$normal[,2], UpperAIC=boot.ci(analysis_boot, type="norm", index=1)$normal[,3], BIC=mean(analysis_boot$t[,2]), LowerBIC=boot.ci(analysis_boot, type="norm", index=2)$normal[,2], UpperBIC=boot.ci(analysis_boot, type="norm", index=2)$normal[,3]))


### MULTIPLE IMPUTATION ###

# #store AIC
# results = data.frame(Model=NULL, AIC=NULL, LowerAIC=NULL, UpperAIC=NULL, BIC=NULL, LowerBIC=NULL, UpperBIC=NULL)
# 
# #imputation
# analysis_mi = mice(analysis, m=5)
# 
# #retrieve each data set and combine into single dataset for bootstrapping
# analysis_mi_data = rbind(complete(analysis_mi, action=1), complete(analysis_mi, action=2), complete(analysis_mi, action=3), complete(analysis_mi, action=4), complete(analysis_mi, action=5))
# 
# #UTD
# 
# #regression modeling
# regressions = with(data=analysis_mi,exp=glm(case~as.factor(pertussis_vax_utd)+as.factor(race)+as.factor(mother_foreign_born)+as.factor(mother_insurance),family=binomial(link="logit")))
# 
# #summary
# summary(pool(regressions))
# 
# #boot strap the 95% CI and mean for AIC & BIC
# analysis_mi_boot = boot(analysis_mi_data, bootAIC, R=1000, model=1)
# #mean(analysis_mi_boot$t[,3])
# results = rbind(results, data.frame(Model="UTD", AIC=mean(analysis_mi_boot$t[,1]), LowerAIC=boot.ci(analysis_mi_boot, type="norm", index=1)$normal[,2], UpperAIC=boot.ci(analysis_mi_boot, type="norm", index=1)$normal[,3], BIC=mean(analysis_mi_boot$t[,2]), LowerBIC=boot.ci(analysis_mi_boot, type="norm", index=2)$normal[,2], UpperBIC=boot.ci(analysis_mi_boot, type="norm", index=2)$normal[,3]))
# 
# #Delayed UTD
# 
# #regression modeling
# regressions = with(data=analysis_mi,exp=glm(case~as.factor(pertussis_vax_utd_delayed)+as.factor(race)+as.factor(mother_foreign_born)+as.factor(mother_insurance),family=binomial(link="logit")))
# 
# #summary
# summary(pool(regressions))
# 
# #boot strap the 95% CI and mean for AIC & BIC
# analysis_mi_boot = boot(analysis_mi_data, bootAIC, R=1000, model=2)
# #mean(analysis_mi_boot$t[,3])
# results = rbind(results, data.frame(Model="Delayed UTD", AIC=mean(analysis_mi_boot$t[,1]), LowerAIC=boot.ci(analysis_mi_boot, type="norm", index=1)$normal[,2], UpperAIC=boot.ci(analysis_mi_boot, type="norm", index=1)$normal[,3], BIC=mean(analysis_mi_boot$t[,2]), LowerBIC=boot.ci(analysis_mi_boot, type="norm", index=2)$normal[,2], UpperBIC=boot.ci(analysis_mi_boot, type="norm", index=2)$normal[,3]))
# 
# #Categorical shots
# 
# #regression modeling
# regressions = with(data=analysis_mi,exp=glm(case~as.factor(pertussis_vax)+as.factor(race)+as.factor(mother_foreign_born)+as.factor(mother_insurance),family=binomial(link="logit")))
# 
# #summary
# summary(pool(regressions))
# 
# #boot strap the 95% CI and mean for AIC & BIC
# analysis_mi_boot = boot(analysis_mi_data, bootAIC, R=1000, model=3)
# #mean(analysis_mi_boot$t[,3])
# results = rbind(results, data.frame(Model="Categorical doses", AIC=mean(analysis_mi_boot$t[,1]), LowerAIC=boot.ci(analysis_mi_boot, type="norm", index=1)$normal[,2], UpperAIC=boot.ci(analysis_mi_boot, type="norm", index=1)$normal[,3], BIC=mean(analysis_mi_boot$t[,2]), LowerBIC=boot.ci(analysis_mi_boot, type="norm", index=2)$normal[,2], UpperBIC=boot.ci(analysis_mi_boot, type="norm", index=2)$normal[,3]))
# 
# #Continuous shots
# 
# #regression modeling
# regressions = with(data=analysis_mi,exp=glm(case~pertussis_vax+as.factor(race)+as.factor(mother_foreign_born)+as.factor(mother_insurance),family=binomial(link="logit")))
# 
# #summary
# summary(pool(regressions))
# 
# #boot strap the 95% CI and mean for AIC & BIC
# analysis_mi_boot = boot(analysis_mi_data, bootAIC, R=1000, model=4)
# #mean(analysis_mi_boot$t[,3])
# results = rbind(results, data.frame(Model="Continuous doses", AIC=mean(analysis_mi_boot$t[,1]), LowerAIC=boot.ci(analysis_mi_boot, type="norm", index=1)$normal[,2], UpperAIC=boot.ci(analysis_mi_boot, type="norm", index=1)$normal[,3], BIC=mean(analysis_mi_boot$t[,2]), LowerBIC=boot.ci(analysis_mi_boot, type="norm", index=2)$normal[,2], UpperBIC=boot.ci(analysis_mi_boot, type="norm", index=2)$normal[,3]))
# 
# #UTD + Categorical shots
# 
# #regression modeling
# regressions = with(data=analysis_mi,exp=glm(case~as.factor(pertussis_vax_utd)+as.factor(pertussis_vax)+as.factor(race)+as.factor(mother_foreign_born)+as.factor(mother_insurance),family=binomial(link="logit")))
# 
# #summary
# summary(pool(regressions))
# 
# #boot strap the 95% CI and mean for AIC & BIC
# analysis_mi_boot = boot(analysis_mi_data, bootAIC, R=1000, model=5)
# #mean(analysis_mi_boot$t[,3])
# results = rbind(results, data.frame(Model="UTD +\nCategorical doses", AIC=mean(analysis_mi_boot$t[,1]), LowerAIC=boot.ci(analysis_mi_boot, type="norm", index=1)$normal[,2], UpperAIC=boot.ci(analysis_mi_boot, type="norm", index=1)$normal[,3], BIC=mean(analysis_mi_boot$t[,2]), LowerBIC=boot.ci(analysis_mi_boot, type="norm", index=2)$normal[,2], UpperBIC=boot.ci(analysis_mi_boot, type="norm", index=2)$normal[,3]))
# 
# #Delayed UTD + Categorical shots
# 
# #regression modeling
# regressions = with(data=analysis_mi,exp=glm(case~as.factor(pertussis_vax_utd_delayed)+as.factor(pertussis_vax)+as.factor(race)+as.factor(mother_foreign_born)+as.factor(mother_insurance),family=binomial(link="logit")))
# 
# #summary
# summary(pool(regressions))
# 
# #boot strap the 95% CI and mean for AIC & BIC
# analysis_mi_boot = boot(analysis_mi_data, bootAIC, R=1000, model=6)
# #mean(analysis_mi_boot$t[,3])
# results = rbind(results, data.frame(Model="Delayed UTD +\nCategorical doses", AIC=mean(analysis_mi_boot$t[,1]), LowerAIC=boot.ci(analysis_mi_boot, type="norm", index=1)$normal[,2], UpperAIC=boot.ci(analysis_mi_boot, type="norm", index=1)$normal[,3], BIC=mean(analysis_mi_boot$t[,2]), LowerBIC=boot.ci(analysis_mi_boot, type="norm", index=2)$normal[,2], UpperBIC=boot.ci(analysis_mi_boot, type="norm", index=2)$normal[,3]))


### ANALYSIS: plot results ###

#reverse order for plotting
results = results[nrow(results):1, ]
rownames(results) = NULL
  
#plot AICs
plot(results[, "AIC"], 1:nrow(results), xlim=c(850,1100), ylab="", axes=F, main="a. Vaccination coverage candidate models by AIC", xlab="Mean AIC and 95% CI", pch=18)
axis(side=1, pos=0.6) 
axis(side=2, las=1, at=1:nrow(results), cex.axis=0.6, labels=(results[,"Model"]), pos=850)

#add confidence bands
for (i in 1:nrow(results))
{
  arrows(as.numeric(results[i, "AIC"]), i, as.numeric(results[i, "LowerAIC"]), i, angle=90, length=0.1)
  arrows(as.numeric(results[i, "AIC"]), i, as.numeric(results[i, "UpperAIC"]), i, angle=90, length=0.1)  
}
rm(i)

#relative to UTD
abline(v=as.numeric(results[results$Model=="UTD","AIC"]), lty=2)

#plot BICs
plot(results[, "BIC"], 1:nrow(results), xlim=c(900,1100), ylab="", axes=F, main="b. Vaccination coverage candidate models by BIC", xlab="Mean BIC and 95% CI", pch=18)
axis(side=1, pos=0.6) 
axis(side=2, las=1, at=1:nrow(results), cex.axis=0.6, labels=(results[,"Model"]), pos=900)

#add confidence bands
for (i in 1:nrow(results))
{
  arrows(as.numeric(results[i, "BIC"]), i, as.numeric(results[i, "LowerBIC"]), i, angle=90, length=0.1)
  arrows(as.numeric(results[i, "BIC"]), i, as.numeric(results[i, "UpperBIC"]), i, angle=90, length=0.1)  
}
rm(i)

#relative to UTD
abline(v=as.numeric(results[results$Model=="UTD","BIC"]), lty=2)


### ANALYSIS: PAR ###

#this assumes crude, not adjusted, but minimal confounding in UTD model
prevUTD = sum(analysis$pertussis_vax_utd, na.rm=T) / length(na.omit(analysis$pertussis_vax_utd))
RR = 0.48
(prevUTD*(RR-1))/(prevUTD*(RR-1)+1)*100
