###################
##               ##
##  Vine Copula  ##
##               ##
##               ##
###################


## Libraries

library(tidyr)
library(dplyr)
library(rugarch)
library(ggplot2)
library(VineCopula)
library(rvinecopulib)
library(MASS)
library(rafalib)
library(readr)
library(fGarch)
library(fable)
library(tsibble)
library(tsibbledata)
library(feasts)
library(readr)
library(tidyverse) 
library(forecast)
library(qqplotr)
library(gridExtra)
library(statmod)
library(tweedie)
library(gamlss)
library(gamlss.inf)
library(gamlss.tr)
library(rugarch)
library(MTS)
library(fitdistrplus)
library(dgof)
library(gamlss.dist)
library(gamlss.add)
library(TSP)
library(smooth)
library(forecast)
library(scoringutils)
library(kdecopula)

## Fit the models
##
# Maximum likelihood estimation of the Tweedie index parameter p
out_Go <- tweedie.profile( UK_data_tsi$Google~1, do.plot=TRUE, verbose=TRUE)
# Fit the glm
fit_Go_Tweedie <- glm( Google ~ date, 
                       data=UK_data_tsi,
                       family=tweedie(var.power=out_Go$xi.max, link.power=0) )
#
fit_NA_SH2 <- gamlss( newAdmissions ~ date, 
                      data=UK_data_tsi,
                      family=SHASHo2,
                      control=gamlss.control(n.cyc=1500))
#
fit_OB_SHASHo<-gamlss(covidOccupiedMVBeds~date,
                      data=UK_data_tsi,
                      family = SHASHo,
                      control = gamlss.control(n.cyc = 1500))
#
fit_deaths_SHo <- gamlss( deaths_count ~ date, 
                          data=UK_data_tsi,
                          family=SHASHo,
                          control=gamlss.control(n.cyc=1500))
#
fit_TW_SEP4 <- gamlss( Total_tweets ~ date, 
                       data=UK_data_tsi,
                       family=SEP4,
                       control=gamlss.control(n.cyc=1500))
#
fit_Bi_net <- gamlss( Bing ~ date, 
                      data=UK_data_tsi,
                      family=NET,
                      control = gamlss.control(c.crit = 0.1,
                                               n.cyc = 100))
#
fit_Af_sst <- gamlss( Afinn ~ date, 
                      data=UK_data_tsi,
                      family=SST,
                      control = gamlss.control(c.crit = 0.1,
                                               n.cyc = 100))
#
VTadj<-UK_data_tsi$newVirusTests/1000
UK_data_tsi<-cbind(UK_data_tsi,VTadj)
#
spec_VT <- ugarchspec(mean.model=list(armaOrder=c(2,1), arfima = FALSE),
                      variance.model=list(garchOrder=c(1,1)),
                      distribution.model="std") 
fit_VT_ar_ga <- ugarchfit(spec = spec_VT, data =UK_data_tsi$VTadj,solver="hybrid" )
#
spec_HC <- ugarchspec(mean.model=list(armaOrder=c(1,1), arfima = FALSE),
                      variance.model=list(garchOrder=c(2,1)),
                      distribution.model="std") 
fit_HC_ar_ga <- ugarchfit(spec = spec_HC, data =UK_data_tsi$hospitalCases,solver="hybrid" )
#
spec_NC <- ugarchspec(mean.model=list(armaOrder=c(1,2), arfima = FALSE),
                      variance.model=list(garchOrder=c(1,1)),
                      distribution.model="std") 
fit_NC_ar_ga <- ugarchfit(spec = spec_NC, data =UK_data_tsi$Cases,solver="hybrid" )
##
#
#
## get the residuals
r.Go<-residuals(fit_Go_Tweedie,standardize=TRUE,sd=sd(fit_Go_Tweedie$residuals))
r.NA<-residuals(fit_NA_SH2,standardize=TRUE)
r.OB<-residuals(fit_OB_SHASHo,standardize=TRUE)
r.Deaths<-residuals(fit_deaths_SHo,standardize=TRUE)
r.TW<-residuals(fit_TW_SEP4,standardize=TRUE)
r.Bing<-residuals(fit_Bi_net,standardize=TRUE)
r.Afinn<-residuals(fit_Af_sst,standardize=TRUE)
r.VT<-residuals(fit_VT_ar_ga,standardize=TRUE)
r.HC<-residuals(fit_HC_ar_ga,standardize=TRUE)
r.NC<-residuals(fit_NC_ar_ga,standardize=TRUE)

sd(r.Bing)
sd(r.Afinn)

u.Go<-pnorm(r.Go,sd=sd(r.Go))
u.NA<-pnorm(r.NA,sd=sd(r.NA))
u.OB<-pnorm(r.OB)
u.Deaths<-pnorm(r.Deaths,sd=sd(r.Deaths))
u.TW<-pnorm(r.TW,sd=sd(r.TW))
u.Bing<-pnorm(r.Bing)
u.Afinn<-pnorm(r.Afinn,sd=sd(r.Afinn))
u.VT<-pnorm(r.VT,sd=sd(r.VT))
u.HC<-pstd(r.HC,sd=sd(r.HC))
u.NC<-pstd(r.NC)
u<-cbind(u.Go,u.NA,u.OB,u.Deaths,u.TW,u.Bing,u.Afinn,u.VT,u.HC,u.NC)

u<-as.copuladata(as.matrix(u))
pairs(u)

## R vines

rv1<-vinecop(u,family_set="par",presel=TRUE)

# Contour plots
contour(rv1,cex.nums=.8)
# Looks like a few independence copulas, particularly in the lower trees.

# Fit results
sum.rv1<- summary(rv1)
print(as.data.frame(sum.rv1), digits = 2)

## Predictions with rolling window
#
## tot number of observations
Tot <- nrow(UK_data_tsi)
Tot
#
## training observations - 90%
R <- round(Tot*90/100)
R
#
## test observations - 10%
P <- round(Tot*10/100)
P
#
## fit copula based on rolling window 
## based on interval [t-R,t-1]
#
### transform residuals into pseudo-observations 
### using parametric method
#
sim<-100
#
# create empty matrix to store simulated copula u-data
UK_sim <- array(0, c((P),10,100))
#
#
# Pull out the residuals for each model
# puting in a different format
res_Go <- fit_Go_Tweedie %>% residuals(standard=TRUE) %>% as_tibble() %>% as.data.frame()
res_NA<-fit_NA_SH2 %>% residuals(standard=TRUE) %>% as_tibble() %>% as.data.frame()
res_OB<-fit_OB_SHASHo %>% residuals(standard=TRUE) %>% as_tibble() %>% as.data.frame()
res_Deaths<-fit_deaths_SHo %>% residuals(standard=TRUE) %>% as_tibble() %>% as.data.frame()
res_TW<-fit_TW_SEP4 %>% residuals(standard=TRUE) %>% as_tibble() %>% as.data.frame()
res_Bing<-fit_Bi_net %>% residuals(standard=TRUE) %>% as_tibble() %>% as.data.frame()
res_Afinn<-fit_Af_sst %>% residuals(standard=TRUE) %>% as_tibble() %>% as.data.frame()
res_VT<-fit_VT_ar_ga %>% residuals(standard=TRUE) %>% as_tibble() %>% as.data.frame()
res_HC<-fit_HC_ar_ga %>% residuals(standard=TRUE) %>% as_tibble() %>% as.data.frame()
res_NC<-fit_NC_ar_ga %>% residuals(standard=TRUE) %>% as_tibble() %>% as.data.frame()
#

for (t in 1:P){
  #
  Admissions <- pnorm(res_NA$value[t:(R+t)],sd=sd(res_NA$value[t:(R+t)]))
  Afinn <- pnorm(res_Afinn$value[t:(R+t)],sd=sd(res_Afinn$value[t:(R+t)]))
  Bing <- pnorm(res_Bing$value[t:(R+t)])
  Cases<-pstd(res_NC$V1[t:(R+t)])
  Deaths <- pnorm(res_Deaths$value[t:(R+t)],sd=sd(res_Deaths$value[t:(R+t)]))
  Google <- pnorm(res_Go$value[t:(R+t)],sd=sd(res_Go$value[t:(R+t)]))
  Hospital <- pstd(res_HC$V1[t:(R+t)],sd=sd(res_HC$V1[t:(R+t)]))
  ICU_Beds <- pnorm(res_OB$value[t:(R+t)])
  Tweets <- pnorm(res_TW$value[t:(R+t)],sd=sd(res_TW$value[t:(R+t)]))
  VirusTests <- pnorm(res_VT$V1[t:(R+t)],sd=sd(res_VT$V1[t:(R+t)]))
  #
  #
  # create a data matrix
  udata_pa <- cbind(Admissions,Afinn,Bing,Cases,Deaths,Google,Hospital,ICU_Beds,Tweets,VirusTests)
  #
  # transform into copula data
  udata_pa <- as.copuladata(udata_pa)
  #
  #
  # select the R-vine structure, families and parameters
  # we allow for the copula families: Gauss, t, Clayton, Gumbel, Frank and Joe and their rotations
  RVM_pa <- RVineStructureSelect(udata_pa, familyset=NA)
  #
  # simulate one u-data vector for each time point
  UK_sim[t,,] <- RVineSim(N=sim, RVM=RVM_pa)
  #
  print(t)
}
#
# check results for simulated u-data
head(UK_sim)
#
#
###
### calculate predicted values for each marginal
###
#
### Google trends
#
pred_Go <- predict(fit_Go_Tweedie, type="response")[(R+1):(Tot)]
head(pred_Go)
#
Dates_test <- UK_data_tsi$date[(R+1):(Tot)]
#
data_test <- data.frame(Dates_test,pred_Go)
#
#
### Twitter
#
pred_Tw <- predict(fit_TW_SEP4,type="response")[(R+1):(Tot)]
head(pred_Tw)
#
#
data_test1 <- data.frame(Dates_test,pred_Tw)
#
#
### Deaths
#
pred_deaths<-predict(fit_deaths_SHo,type="response")[(R+1):(Tot)]
head(pred_deaths)
#
#
data_test2<-data.frame(Dates_test,pred_deaths)
#
#
### Bing
#
pred_Bing<-predict(fit_Bi_net,type="response")[(R+1):(Tot)]
head(pred_Bing)
#
#
data_test3<-data.frame(Dates_test,pred_Bing)
#
#
### Afinn
pred_Afinn<-predict(fit_Af_sst,type="response")[(R+1):(Tot)]
head(pred_Afinn)
#
#
data_test4<-data.frame(Dates_test,pred_Afinn)
#
#
### New Admissions
pred_NA<-predict(fit_NA_SH2,type="response")[(R+1):(Tot)]
head(pred_NA)
#
#
data_test5<-data.frame(Dates_test,pred_NA)
#
#
### Occupied Beds
#
#
pred_OB<-predict(fit_OB_SHASHo,type="response")[(R+1):(Tot)]
head(pred_OB)
#
#
data_test6<-data.frame(Dates_test,pred_OB)
#
#
### Hospital Cases
#
HC_pred<-ugarchforecast(fit_HC_ar_ga,n.ahead=38,data = UK_data_ren[1:346,5,drop=FALSE],out.sample = 38) 
pred2<-fitted(HC_pred) # Has an option of n.roll but it keeps throwing me an error
pred_HC<-pred2[,1]

data_test7 <- data.frame(Dates_test,pred_HC)
#
#
### Cases
#
NC_pred<-ugarchforecast(fit_NC_ar_ga,n.ahead = 38,data = UK_data_ren[1:346,7,drop=FALSE],out.sample = 38)
pred3<-fitted(NC_pred)
pred_NC<-pred3[,1]
#
#
data_test8<-data.frame(Dates_test,pred_NC)
#
#
### Virus Tests per 1000
#
#
UK_data_ren<-cbind(UK_data_ren,VTadj)
VT_pred<-ugarchforecast(fit_VT_ar_ga,n.ahead = 38,data = UK_data_ren[1:346,12,drop=FALSE],out.sample = 38)
pred4<-fitted(VT_pred)
pred_VT<-pred4[,1]
#
#
data_test9<-data.frame(Dates_test,pred_VT)
#
#

Vector<-c(data_test$pred_Go[1:37],data_test5$pred_NA[1:37],data_test6$pred_OB[1:37],
          data_test2$pred_deaths[1:37],data_test1$pred_Tw[1:37],data_test3$pred_Bing[1:37],
          data_test4$pred_Afinn[1:37],data_test9$pred_VT[1:37],data_test7$pred_HC[1:37],
          data_test8$pred_NC[1:37])

pred_data<-array(c(Vector),c(37,10,100))

pred_data

results<-qnorm(UK_sim) +pred_data

# Take the average of all the simulations
UK_sim_avg <- apply(results, c(1,2), mean)

# 
#
Dates_test <- UK_data_tsi$date[(R+1):(Tot-1)]
#
UK_sim_avg<-data.frame(Dates_test,UK_sim_avg)
names(UK_sim_avg)[2] <- "Google"
names(UK_sim_avg)[3] <- "New_Admissions"
names(UK_sim_avg)[4] <- "Occupied_Beds"
names(UK_sim_avg)[5] <- "Deaths"
names(UK_sim_avg)[6] <- "Twitter"
names(UK_sim_avg)[7] <- "Bing"
names(UK_sim_avg)[8] <- "Afinn"
names(UK_sim_avg)[9] <- "Virus_Tests"
names(UK_sim_avg)[10] <- "Hospital_cases"
names(UK_sim_avg)[11] <- "New_Cases"

#
#
#
## Predictions using just a Gaussian copula

UK_sim_G <- array(0, c((P-1),10,sim))

for (t in 1:P){
  #
  UK_1_pa <- pnorm(res_Go$value[t:(R+t)],sd=sd(res_Go$value[t:(R+t)]))
  UK_2_pa <- pnorm(res_NA$value[t:(R+t)],sd=sd(res_NA$value[t:(R+t)]))
  UK_3_pa <- pnorm(res_OB$value[t:(R+t)])
  UK_4_pa <- pnorm(res_Deaths$value[t:(R+t)],sd=sd(res_Deaths$value[t:(R+t)]))
  UK_5_pa <- pnorm(res_TW$value[t:(R+t)],sd=sd(res_TW$value[t:(R+t)]))
  UK_6_pa <- pnorm(res_Bing$value[t:(R+t)])
  UK_7_pa <- pnorm(res_Afinn$value[t:(R+t)],sd=sd(res_Afinn$value[t:(R+t)]))
  UK_8_pa <- pnorm(res_VT$V1[t:(R+t)],sd=sd(res_VT$V1[t:(R+t)]))
  UK_9_pa <- pstd(res_HC$V1[t:(R+t)],sd=sd(res_HC$V1[t:(R+t)]))
  UK_10_pa<-pstd(res_NA$value[t:(R+t)])
  #
  #
  # create a data matrix
  udata_pa <- cbind(UK_1_pa,UK_2_pa,UK_3_pa,UK_4_pa,UK_5_pa,
                    UK_6_pa,UK_7_pa,UK_8_pa,UK_9_pa,UK_10_pa)
  #
  # transform into copula data
  udata_pa <- as.copuladata(udata_pa)
  #
  #
  # select the R-vine structure, families and parameters
  # we allow for the copula families: Gauss, t, Clayton, Gumbel, Frank and Joe and their rotations
  RVM_pa_G <- RVineStructureSelect(udata_pa, familyset=1)
  #
  # simulate one u-data vector for each time point
  UK_sim_G[t,,] <- RVineSim(N=sim, RVM=RVM_pa_G)
  #
  print(t)
}
#
# check results for simulated u-data
head(UK_sim_G)

#
results_G<-qnorm(UK_sim_G) +pred_data
# 
# Take the average of all the simulations
UK_sim_G_avg <- apply(results_G, c(1,2), mean)
#
# Take the average of all the simulations
UK_sim_G_avg<-data.frame(Dates_test,UK_sim_G_avg)
names(UK_sim_G_avg)[2] <- "Google"
names(UK_sim_G_avg)[3] <- "New_Admissions"
names(UK_sim_G_avg)[4] <- "Occupied_Beds"
names(UK_sim_G_avg)[5] <- "Deaths"
names(UK_sim_G_avg)[6] <- "Twitter"
names(UK_sim_G_avg)[7] <- "Bing"
names(UK_sim_G_avg)[8] <- "Afinn"
names(UK_sim_G_avg)[9] <- "Virus_Tests"
names(UK_sim_G_avg)[10] <- "Hospital_cases"
names(UK_sim_G_avg)[11] <- "New_Cases"
# 
#
###
### calculate predicted values for each marginal
###
#

## Predictions using just the independence copula


UK_sim_I <- array(0, c((P),10,sim))

for (t in 1:P){
  #
  UK_1_pa <- pnorm(res_Go$value[t:(R+t)],sd=sd(res_Go$value[t:(R+t)]))
  UK_2_pa <- pnorm(res_NA$value[t:(R+t)],sd=sd(res_NA$value[t:(R+t)]))
  UK_3_pa <- pnorm(res_OB$value[t:(R+t)])
  UK_4_pa <- pnorm(res_Deaths$value[t:(R+t)],sd=sd(res_Deaths$value[t:(R+t)]))
  UK_5_pa <- pnorm(res_TW$value[t:(R+t)],sd=sd(res_TW$value[t:(R+t)]))
  UK_6_pa <- pnorm(res_Bing$value[t:(R+t)])
  UK_7_pa <- pnorm(res_Afinn$value[t:(R+t)],sd=sd(res_Afinn$value[t:(R+t)]))
  UK_8_pa <- pnorm(res_VT$V1[t:(R+t)],sd=sd(res_VT$V1[t:(R+t)]))
  UK_9_pa <- pstd(res_HC$V1[t:(R+t)],sd=sd(res_HC$V1[t:(R+t)]))
  UK_10_pa<-pstd(res_NA$value[t:(R+t)])
  #
  #
  # create a data matrix
  udata_pa <- cbind(UK_1_pa,UK_2_pa,UK_3_pa,UK_4_pa,UK_5_pa,
                    UK_6_pa,UK_7_pa,UK_8_pa,UK_9_pa,UK_10_pa)
  #
  # transform into copula data
  udata_pa <- as.copuladata(udata_pa)
  #
  #
  # select the R-vine structure, families and parameters
  # we allow for the copula families: Gauss, t, Clayton, Gumbel, Frank and Joe and their rotations
  RVM_pa_I <- RVineStructureSelect(udata_pa, familyset=0)
  #
  # simulate one u-data vector for each time point
  UK_sim_I[t,,] <- RVineSim(N=sim, RVM=RVM_pa_I)
  #
  print(t)
}
#
# check results for simulated u-data
head(UK_sim_I)
#
results_I<-qnorm(UK_sim_I) +pred_data
# 
# Take the average of all the simulations
UK_sim_I_avg <- apply(results_I, c(1,2), mean)
#
# Take the average of all the simulations
UK_sim_I_avg<-data.frame(Dates_test,UK_sim_I_avg)
names(UK_sim_I_avg)[2] <- "Google"
names(UK_sim_I_avg)[3] <- "New_Admissions"
names(UK_sim_I_avg)[4] <- "Occupied_Beds"
names(UK_sim_I_avg)[5] <- "Deaths"
names(UK_sim_I_avg)[6] <- "Twitter"
names(UK_sim_I_avg)[7] <- "Bing"
names(UK_sim_I_avg)[8] <- "Afinn"
names(UK_sim_I_avg)[9] <- "Virus_Tests"
names(UK_sim_I_avg)[10] <- "Hospital_cases"
names(UK_sim_I_avg)[11] <- "New_Cases"
# 
#

## Google
test_Go<-UK_data_ren[347:384,c(1,2)]

## Twitter
test_Tw<-UK_data_ren[347:384,c(1,11)]

## Bing
test_Bing<-UK_data_ren[347:384,c(1,9)]

## Afinn
test_Afinn<-UK_data_ren[347:384,c(1,10)]

## Virus Test
test_VT<-UK_data_ren[347:384,c(1,12)]

## New admissions
test_NA<-UK_data_ren[347:384,c(1,4)]

## Hospital Cases
test_HC<-UK_data_ren[347:384,c(1,5)]

## Occupied Beds
test_OB<-UK_data_ren[347:384,c(1,6)]

## Cases
test_NC<-UK_data_ren[347:384,c(1,7)]

## Deaths
test_deaths<-UK_data_ren[347:384,c(1,8)]

### MSE

## Google
MSE(test_Go$Google[1:37],UK_sim_avg$Google)
MSE(test_Go$Google[1:37],UK_sim_G_avg$Google)
MSE(test_Go$Google[1:37],UK_sim_I_avg$Google)

## Twitter
MSE(test_Tw$Total_tweets[1:37],UK_sim_avg$Twitter)
MSE(test_Tw$Total_tweets[1:37],UK_sim_G_avg$Twitter)
MSE(test_Tw$Total_tweets[1:37],UK_sim_I_avg$Twitter)

## Bing
MSE(test_Bing$Bing[1:37],UK_sim_avg$Bing)
MSE(test_Bing$Bing[1:37],UK_sim_G_avg$Bing)
MSE(test_Bing$Bing[1:37],UK_sim_I_avg$Bing)

## Afinn
MSE(test_Afinn$Afinn[1:37],UK_sim_avg$Afinn)
MSE(test_Afinn$Afinn[1:37],UK_sim_G_avg$Afinn)
MSE(test_Afinn$Afinn[1:37],UK_sim_I_avg$Afinn)

## Virus tests
MSE(test_VT$VTadj[1:37],UK_sim_avg$Virus_Tests)
MSE(test_VT$VTadj[1:37],UK_sim_G_avg$Virus_Tests)
MSE(test_VT$VTadj[1:37],UK_sim_I_avg$Virus_Tests)

## New admissions
MSE(test_NA$newAdmissions[1:37],UK_sim_avg$New_Admissions)
MSE(test_NA$newAdmissions[1:37],UK_sim_G_avg$New_Admissions)
MSE(test_NA$newAdmissions[1:37],UK_sim_I_avg$New_Admissions)

## Hospital cases
MSE(test_HC$hospitalCases[1:37],UK_sim_avg$Hospital_cases)
MSE(test_HC$hospitalCases[1:37],UK_sim_G_avg$Hospital_cases)
MSE(test_HC$hospitalCases[1:37],UK_sim_I_avg$Hospital_cases)

## Occupied beds
MSE(test_OB$covidOccupiedMVBeds[1:37],UK_sim_avg$Occupied_Beds)
MSE(test_OB$covidOccupiedMVBeds[1:37],UK_sim_G_avg$Occupied_Beds)
MSE(test_OB$covidOccupiedMVBeds[1:37],UK_sim_I_avg$Occupied_Beds)

## Cases
MSE(test_NC$Cases[1:37],UK_sim_avg$New_Cases)
MSE(test_NC$Cases[1:37],UK_sim_G_avg$New_Cases)
MSE(test_NC$Cases[1:37],UK_sim_I_avg$New_Cases)

## Deaths
MSE(test_deaths$deaths_count[1:37],UK_sim_avg$Deaths)
MSE(test_deaths$deaths_count[1:37],UK_sim_G_avg$Deaths)
MSE(test_deaths$deaths_count[1:37],UK_sim_I_avg$Deaths)

### Interval score

interval_range = 90

## Google

x_bar_Go<-mean(test_Go$Google)
sd_Go<-sd(test_Go$Google)
Go_res_M<-test_Go$Google[1:37] - UK_sim_avg$Google
res_sum_sq_M_Go<-sum(Go_res_M^2)
sigma_Go_M<-sqrt((1/(37-2-1))*res_sum_sq_M_Go)
value_M_Go<-1.96*sigma_Go_M*sqrt(1+(1/37)+(UK_sim_avg$Google-x_bar_Go)^2/((37-1)*sd_Go^2))

lower_M_Go<-UK_sim_avg$Google-value_M_Go
upper_M_Go<-UK_sim_avg$Google+value_M_Go

mean(interval_score(true_values = test_Go$Google[1:37],
                    lower = lower_M_Go,
                    upper = upper_M_Go,
                    interval_range = interval_range))

Go_res_G<-test_Go$Google[1:37] - UK_sim_G_avg$Google
res_sum_sq_G_Go<-sum(Go_res_G^2)
sigma_Go_G<-sqrt((1/(37-2-1))*res_sum_sq_G_Go)
value_G_Go<-1.96*sigma_Go_G*sqrt(1+(1/37)+(UK_sim_G_avg$Google-x_bar_Go)^2/((37-1)*sd_Go^2))


lower_G_GO<-UK_sim_G_avg$Google-value_G_Go
upper_G_GO<-UK_sim_G_avg$Google+value_G_Go

mean(interval_score(true_values = test_Go$Google[1:37],
                    lower = lower_G_GO,
                    upper = upper_G_GO,
                    interval_range = interval_range))

Go_res_I<-test_Go$Google[1:37] - UK_sim_I_avg$Google
res_sum_sq_I_Go<-sum(Go_res_I^2)
sigma_Go_I<-sqrt((1/(37-2-1))*res_sum_sq_I_Go)
value_I_Go<-1.96*sigma_Go_I*sqrt(1+(1/37)+(UK_sim_I_avg$Google-x_bar_Go)^2/((37-1)*sd_Go^2))

lower_I_GO<-UK_sim_I_avg$Google-value_I_Go
upper_I_GO<-UK_sim_I_avg$Google+value_I_Go

mean(interval_score(true_values = test_Go$Google[1:37],
                    lower = lower_I_GO,
                    upper = upper_I_GO,
                    interval_range = interval_range))

## Twitter

x_bar_Tw<-mean(test_Tw$Total_tweets)
sd_Tw<-sd(test_Tw$Total_tweets)
Tw_res_M<-test_Tw$Total_tweets[1:37] - UK_sim_avg$Twitter
res_sum_sq_M_Tw<-sum(Tw_res_M^2)
sigma_Tw_M<-sqrt((1/(37-4-1))*res_sum_sq_M_Tw)
value_M_Tw<-1.96*sigma_Tw_M*sqrt(1+(1/37)+(UK_sim_avg$Twitter-x_bar_Tw)^2/((37-1)*sd_Tw^2))

lower_M_TW<-UK_sim_avg$Twitter-value_M_Tw
upper_M_TW<-UK_sim_avg$Twitter+value_M_Tw

mean(interval_score(true_values = test_Tw$Total_tweets[1:37],
                    lower = lower_M_TW,
                    upper = upper_M_TW,
                    interval_range = interval_range))

Tw_res_G<-test_Tw$Total_tweets[1:37] - UK_sim_G_avg$Twitter
res_sum_sq_G_Tw<-sum(Tw_res_G^2)
sigma_Tw_G<-sqrt((1/(37-4-1))*res_sum_sq_G_Tw)
value_G_Tw<-1.96*sigma_Tw_G*sqrt(1+(1/37)+(UK_sim_G_avg$Twitter-x_bar_Tw)^2/((37-1)*sd_Tw^2))


lower_G_TW<-UK_sim_G_avg$Twitter-value_G_Tw
upper_G_TW<-UK_sim_G_avg$Twitter+value_G_Tw

mean(interval_score(true_values = test_Tw$Total_tweets[1:37],
                    lower = lower_G_TW,
                    upper = upper_G_TW,
                    interval_range = interval_range))

Tw_res_I<-test_Tw$Total_tweets[1:37] - UK_sim_I_avg$Twitter
res_sum_sq_I_Tw<-sum(Tw_res_I^2)
sigma_Tw_I<-sqrt((1/(37-4-1))*res_sum_sq_I_Tw)
value_I_Tw<-1.96*sigma_Tw_I*sqrt(1+(1/37)+(UK_sim_I_avg$Twitter-x_bar_Tw)^2/((37-1)*sd_Tw^2))

lower_I_TW<-UK_sim_I_avg$Twitter-value_I_Tw
upper_I_TW<-UK_sim_I_avg$Twitter+value_I_Tw

mean(interval_score(true_values = test_Tw$Total_tweets[1:37],
                    lower = lower_I_TW,
                    upper = upper_I_TW,
                    interval_range = interval_range))

## New cases

x_bar_NC<-mean(test_NC$Cases)
sd_NC<-sd(test_NC$Cases)
NC_res_M<-test_NC$Cases[1:37] - UK_sim_avg$New_Cases
res_sum_sq_M_NC<-sum(NC_res_M^2)
sigma_NC_M<-sqrt((1/(37-5-1))*res_sum_sq_M_NC)
value_M_NC<-1.96*sigma_NC_M*sqrt(1+(1/37)+(UK_sim_avg$New_Cases-x_bar_NC)^2/((37-1)*sd_NC^2))


lower_M_NC<-UK_sim_avg$New_Cases-value_M_NC
upper_M_NC<-UK_sim_avg$New_Cases+value_M_NC

mean(interval_score(true_values = test_NC$Cases[1:37],
                    lower = lower_M_NC,
                    upper = upper_M_NC,
                    interval_range = interval_range))

NC_res_G<-test_NC$Cases[1:37] - UK_sim_G_avg$New_Cases
res_sum_sq_G_NC<-sum(NC_res_G^2)
sigma_NC_G<-sqrt((1/(37-5-1))*res_sum_sq_G_NC)
value_G_NC<-1.96*sigma_NC_G*sqrt(1+(1/37)+(UK_sim_G_avg$New_Cases-x_bar_NC)^2/((37-1)*sd_NC^2))

lower_G_NC<-UK_sim_G_avg$New_Cases-value_G_NC
upper_G_NC<-UK_sim_G_avg$New_Cases+value_G_NC

mean(interval_score(true_values = test_NC$Cases[1:37],
                    lower = lower_G_NC,
                    upper = upper_G_NC,
                    interval_range = interval_range))

NC_res_I<-test_NC$Cases[1:37] - UK_sim_I_avg$New_Cases
res_sum_sq_I_NC<-sum(NC_res_I^2)
sigma_NC_I<-sqrt((1/(37-5-1))*res_sum_sq_I_NC)
value_I_NC<-1.96*sigma_NC_I*sqrt(1+(1/37)+(UK_sim_I_avg$New_Cases-x_bar_NC)^2/((37-1)*sd_NC^2))


lower_I_NC<-UK_sim_I_avg$New_Cases-value_I_NC
upper_I_NC<-UK_sim_I_avg$New_Cases+value_I_NC

mean(interval_score(true_values = test_NC$Cases[1:37],
                    lower = lower_I_NC,
                    upper = upper_I_NC,
                    interval_range = interval_range))

## Occupied beds

x_bar_OB<-mean(test_OB$covidOccupiedMVBeds)
sd_OB<-sd(test_OB$covidOccupiedMVBeds)
OB_res_M<-test_OB$covidOccupiedMVBeds[1:37] - UK_sim_avg$Occupied_Beds
res_sum_sq_M<-sum(OB_res_M^2)
sigma_OB_M<-sqrt((1/(37-4-1))*res_sum_sq_M)
value_M<-1.96*sigma_OB_M*sqrt(1+(1/37)+(UK_sim_avg$Occupied_Beds-x_bar_OB)^2/((37-1)*sd_OB^2))

lower_M_OB<-UK_sim_avg$Occupied_Beds-value_M
upper_M_OB<-UK_sim_avg$Occupied_Beds+value_M

mean(interval_score(true_values = test_OB$covidOccupiedMVBeds[1:37],
                    lower = lower_M_OB,
                    upper = upper_M_OB,
                    interval_range = interval_range))

OB_res_G<-test_OB$covidOccupiedMVBeds[1:37] - UK_sim_G_avg$Occupied_Beds
res_sum_sq_G<-sum(OB_res_G^2)
sigma_OB_G<-sqrt((1/(37-4-1))*res_sum_sq_G)
value_G<-1.96*sigma_OB_G*sqrt(1+(1/37)+(UK_sim_G_avg$Occupied_Beds-x_bar_OB)^2/((37-1)*sd_OB^2))

lower_G_OB<-UK_sim_G_avg$Occupied_Beds-value_G
upper_G_OB<-UK_sim_G_avg$Occupied_Beds+value_G

mean(interval_score(true_values = test_OB$covidOccupiedMVBeds[1:37],
                    lower = lower_G_OB,
                    upper = upper_G_OB,
                    interval_range = interval_range))

OB_res_I<-test_OB$covidOccupiedMVBeds[1:37] - UK_sim_I_avg$Occupied_Beds
res_sum_sq_I<-sum(OB_res_I^2)
sigma_OB_I<-sqrt((1/(37-4-1))*res_sum_sq_I)
value_I<-1.96*sigma_OB_I*sqrt(1+(1/37)+(UK_sim_G_avg$Occupied_Beds-x_bar_OB)^2/((37-1)*sd_OB^2))

lower_I_OB<-UK_sim_I_avg$Occupied_Beds-value_I
upper_I_OB<-UK_sim_I_avg$Occupied_Beds+value_I

mean(interval_score(true_values = test_OB$covidOccupiedMVBeds[1:37],
                    lower = lower_I_OB,
                    upper = upper_I_OB,
                    interval_range = interval_range))

## Bing

x_bar_Bi<-mean(test_Bing$Bing)
sd_Bi<-sd(test_Bing$Bing)
Bi_res_M<-test_Bing$Bing[1:37] - UK_sim_avg$Bing
res_sum_sq_M_Bi<-sum(Bi_res_M^2)
sigma_Bi_M<-sqrt((1/(37-4-1))*res_sum_sq_M_Bi)
value_M_Bi<-1.96*sigma_Bi_M*sqrt(1+(1/37)+(UK_sim_avg$Bing-x_bar_Bi)^2/((37-1)*sd_Bi^2))

lower_M_Bi<-UK_sim_avg$Bing-value_M_Bi
upper_M_Bi<-UK_sim_avg$Bing+value_M_Bi

mean(interval_score(true_values = test_Bing$Bing[1:37],
                    lower = lower_M_Bi,
                    upper = upper_M_Bi,
                    interval_range = interval_range))

Bi_res_G<-test_Bing$Bing[1:37] - UK_sim_G_avg$Bing
res_sum_sq_G_Bi<-sum(Bi_res_G^2)
sigma_Bi_G<-sqrt((1/(37-4-1))*res_sum_sq_G_Bi)
value_G_Bi<-1.96*sigma_Bi_G*sqrt(1+(1/37)+(UK_sim_G_avg$Bing-x_bar_Bi)^2/((37-1)*sd_Bi^2))

lower_G_Bi<-UK_sim_G_avg$Bing-value_G_Bi
upper_G_Bi<-UK_sim_G_avg$Bing+value_G_Bi

mean(interval_score(true_values = test_Bing$Bing[1:37],
                    lower = lower_G_Bi,
                    upper = upper_G_Bi,
                    interval_range = interval_range))

Bi_res_I<-test_Bing$Bing[1:37] - UK_sim_I_avg$Bing
res_sum_sq_I_Bi<-sum(Bi_res_I^2)
sigma_Bi_I<-sqrt((1/(37-4-1))*res_sum_sq_I_Bi)
value_I_Bi<-1.96*sigma_Bi_I*sqrt(1+(1/37)+(UK_sim_I_avg$Bing-x_bar_Bi)^2/((37-1)*sd_Bi^2))

lower_I_Bi<-UK_sim_I_avg$Bing-value_I_Bi
upper_I_Bi<-UK_sim_I_avg$Bing+value_I_Bi

mean(interval_score(true_values = test_Bing$Bing[1:37],
                    lower = lower_I_Bi,
                    upper = upper_I_Bi,
                    interval_range = interval_range))
## Afinn

x_bar_Af<-mean(test_Afinn$Afinn)
sd_Af<-sd(test_Afinn$Afinn)
Af_res_M<-test_Afinn$Afinn[1:37] - UK_sim_avg$Afinn
res_sum_sq_M_Af<-sum(Af_res_M^2)
sigma_Af_M<-sqrt((1/(37-4-1))*res_sum_sq_M_Af)
value_M_Af<-1.96*sigma_Af_M*sqrt(1+(1/37)+(UK_sim_avg$Afinn-x_bar_Af)^2/((37-1)*sd_Af^2))

lower_M_Af<-UK_sim_avg$Afinn-value_M_Af
upper_M_Af<-UK_sim_avg$Afinn+value_M_Af

mean(interval_score(true_values = test_Afinn$Afinn[1:37],
                    lower = lower_M_Af,
                    upper = upper_M_Af,
                    interval_range = interval_range))

Af_res_G<-test_Afinn$Afinn[1:37] - UK_sim_G_avg$Afinn
res_sum_sq_G_Af<-sum(Af_res_G^2)
sigma_Af_G<-sqrt((1/(37-4-1))*res_sum_sq_G_Af)
value_G_Af<-1.96*sigma_Af_G*sqrt(1+(1/37)+(UK_sim_G_avg$Afinn-x_bar_Af)^2/((37-1)*sd_Af^2))

lower_G_Af<-UK_sim_G_avg$Afinn-value_G_Af
upper_G_Af<-UK_sim_G_avg$Afinn+value_G_Af

mean(interval_score(true_values = test_Afinn$Afinn[1:37],
                    lower = lower_G_Af,
                    upper = upper_G_Af,
                    interval_range = interval_range))

Af_res_I<-test_Afinn$Afinn[1:37] - UK_sim_I_avg$Afinn
res_sum_sq_I_Af<-sum(Af_res_I^2)
sigma_Af_I<-sqrt((1/(37-4-1))*res_sum_sq_I_Af)
value_I_Af<-1.96*sigma_Af_I*sqrt(1+(1/37)+(UK_sim_I_avg$Afinn-x_bar_Af)^2/((37-1)*sd_Af^2))

lower_I_Af<-UK_sim_I_avg$Afinn-value_I_Af
upper_I_Af<-UK_sim_I_avg$Afinn+value_I_Af

mean(interval_score(true_values = test_Afinn$Afinn[1:37],
                    lower = lower_I_Af,
                    upper = upper_I_Af,
                    interval_range = interval_range))

## Deaths

x_bar_Deaths<-mean(test_deaths$deaths_count)
sd_Deaths<-sd(test_deaths$deaths_count)
Deaths_res_M<-test_deaths$deaths_count[1:37] - UK_sim_avg$Deaths
res_sum_sq_M_Deaths<-sum(Deaths_res_M^2)
sigma_Deaths_M<-sqrt((1/(37-4-1))*res_sum_sq_M_Deaths)
value_M_Deaths<-1.96*sigma_Deaths_M*sqrt(1+(1/37)+(UK_sim_avg$Deaths-x_bar_Deaths)^2/((37-1)*sd_Deaths^2))

lower_M_Deaths<-UK_sim_avg$Deaths-value_M_Deaths
upper_M_Deaths<-UK_sim_avg$Deaths+value_M_Deaths

mean(interval_score(true_values = test_deaths$deaths_count[1:37],
                    lower = lower_M_Deaths,
                    upper = upper_M_Deaths,
                    interval_range = interval_range))

Deaths_res_G<-test_deaths$deaths_count[1:37] - UK_sim_G_avg$Deaths
res_sum_sq_G_Deaths<-sum(Deaths_res_G^2)
sigma_Deaths_G<-sqrt((1/(37-4-1))*res_sum_sq_G_Deaths)
value_G_Deaths<-1.96*sigma_Deaths_G*sqrt(1+(1/37)+(UK_sim_G_avg$Deaths-x_bar_Deaths)^2/((37-1)*sd_Deaths^2))


lower_G_Deaths<-UK_sim_G_avg$Deaths-value_G_Deaths
upper_G_Deaths<-UK_sim_G_avg$Deaths+value_G_Deaths

mean(interval_score(true_values = test_deaths$deaths_count[1:37],
                    lower = lower_G_Deaths,
                    upper = upper_G_Deaths,
                    interval_range = interval_range))

Deaths_res_I<-test_deaths$deaths_count[1:37] - UK_sim_I_avg$Deaths
res_sum_sq_I_Deaths<-sum(Deaths_res_I^2)
sigma_Deaths_I<-sqrt((1/(37-4-1))*res_sum_sq_I_Deaths)
value_I_Deaths<-1.96*sigma_Deaths_I*sqrt(1+(1/37)+(UK_sim_I_avg$Deaths-x_bar_Deaths)^2/((37-1)*sd_Deaths^2))

lower_I_Deaths<-UK_sim_I_avg$Deaths-value_I_Deaths
upper_I_Deaths<-UK_sim_I_avg$Deaths+value_I_Deaths

mean(interval_score(true_values = test_deaths$deaths_count[1:37],
                    lower = lower_I_Deaths,
                    upper = upper_I_Deaths,
                    interval_range = interval_range))

## Hospital Cases

x_bar_HC<-mean(test_HC$hospitalCases)
sd_HC<-sd(test_HC$hospitalCases)
HC_res_M<-test_HC$hospitalCases[1:37] - UK_sim_avg$Hospital_cases
res_sum_sq_M_HC<-sum(HC_res_M^2)
sigma_HC_M<-sqrt((1/(37-5-1))*res_sum_sq_M_HC)
value_M_HC<-1.96*sigma_HC_M*sqrt(1+(1/37)+(UK_sim_avg$Hospital_cases-x_bar_HC)^2/((37-1)*sd_HC^2))

lower_M_HC<-UK_sim_avg$Hospital_cases-value_M_HC
upper_M_HC<-UK_sim_avg$Hospital_cases+value_M_HC

mean(interval_score(true_values = test_HC$hospitalCases[1:37],
                    lower = lower_M_HC,
                    upper = upper_M_HC,
                    interval_range = interval_range))

HC_res_G<-test_HC$hospitalCases[1:37] - UK_sim_G_avg$Hospital_cases
res_sum_sq_G_HC<-sum(HC_res_G^2)
sigma_HC_G<-sqrt((1/(37-5-1))*res_sum_sq_G_HC)
value_G_HC<-1.96*sigma_HC_G*sqrt(1+(1/37)+(UK_sim_G_avg$Hospital_cases-x_bar_HC)^2/((37-1)*sd_HC^2))


lower_G_HC<-UK_sim_G_avg$Hospital_cases-value_G_HC
upper_G_HC<-UK_sim_G_avg$Hospital_cases+value_G_HC

mean(interval_score(true_values = test_HC$hospitalCases[1:37],
                    lower = lower_G_HC,
                    upper = upper_G_HC,
                    interval_range = interval_range))

HC_res_I<-test_HC$hospitalCases[1:37] - UK_sim_I_avg$Hospital_cases
res_sum_sq_I_HC<-sum(HC_res_I^2)
sigma_HC_I<-sqrt((1/(37-5-1))*res_sum_sq_I_HC)
value_I_HC<-1.96*sigma_HC_I*sqrt(1+(1/37)+(UK_sim_I_avg$Hospital_cases-x_bar_HC)^2/((37-1)*sd_HC^2))

lower_I_HC<-UK_sim_I_avg$Hospital_cases-value_I_HC
upper_I_HC<-UK_sim_I_avg$Hospital_cases+value_I_HC

mean(interval_score(true_values = test_HC$hospitalCases[1:37],
                    lower = lower_I_HC,
                    upper = upper_I_HC,
                    interval_range = interval_range))

## New admissions

x_bar_NA<-mean(test_NA$newAdmissions)
sd_NA<-sd(test_NA$newAdmissions)
NA_res_M<-test_NA$newAdmissions[1:37] - UK_sim_avg$New_Admissions
res_sum_sq_M_NA<-sum(NA_res_M^2)
sigma_NA_M<-sqrt((1/(37-4-1))*res_sum_sq_M_NA)
value_M_NA<-1.96*sigma_NA_M*sqrt(1+(1/37)+(UK_sim_avg$New_Admissions-x_bar_NA)^2/((37-1)*sd_NA^2))

lower_M_NA<-UK_sim_avg$New_Admissions-value_M_NA
upper_M_NA<-UK_sim_avg$New_Admissions+value_M_NA

mean(interval_score(true_values = test_NA$newAdmissions[1:37],
                    lower = lower_M_NA,
                    upper = upper_M_NA,
                    interval_range = interval_range))

NA_res_G<-test_NA$newAdmissions[1:37] - UK_sim_G_avg$New_Admissions
res_sum_sq_G_NA<-sum(NA_res_G^2)
sigma_NA_G<-sqrt((1/(37-4-1))*res_sum_sq_G_NA)
value_G_NA<-1.96*sigma_NA_G*sqrt(1+(1/37)+(UK_sim_G_avg$New_Admissions-x_bar_NA)^2/((37-1)*sd_NA^2))


lower_G_NA<-UK_sim_G_avg$New_Admissions-value_G_NA
upper_G_NA<-UK_sim_G_avg$New_Admissions+value_G_NA

mean(interval_score(true_values = test_NA$newAdmissions[1:37],
                    lower = lower_G_NA,
                    upper = upper_G_NA,
                    interval_range = interval_range))

NA_res_I<-test_NA$newAdmissions[1:37] - UK_sim_I_avg$New_Admissions
res_sum_sq_I_NA<-sum(NA_res_I^2)
sigma_NA_I<-sqrt((1/(37-4-1))*res_sum_sq_I_NA)
value_I_NA<-1.96*sigma_NA_I*sqrt(1+(1/37)+(UK_sim_I_avg$New_Admissions-x_bar_NA)^2/((37-1)*sd_NA^2))

lower_I_NA<-UK_sim_I_avg$New_Admissions-value_I_NA
upper_I_NA<-UK_sim_I_avg$New_Admissions+value_I_NA

mean(interval_score(true_values = test_NA$newAdmissions[1:37],
                    lower = lower_I_NA,
                    upper = upper_I_NA,
                    interval_range = interval_range))

## Virus tests (adjusted)

x_bar_VT<-mean(test_VT$VTadj)
sd_VT<-sd(test_VT$VTadj)
VT_res_M<-test_VT$VTadj[1:37] - UK_sim_avg$Virus_Tests
res_sum_sq_M_VT<-sum(VT_res_M^2)
sigma_VT_M<-sqrt((1/(37-5-1))*res_sum_sq_M_VT)
value_M_VT<-1.96*sigma_VT_M*sqrt(1+(1/37)+(UK_sim_avg$Virus_Tests-x_bar_VT)^2/((37-1)*sd_VT^2))

lower_M_VT<-UK_sim_avg$Virus_Tests-value_M_VT
upper_M_VT<-UK_sim_avg$Virus_Tests+value_M_VT

mean(interval_score(true_values = test_VT$VTadj[1:37],
                    lower = lower_M_VT,
                    upper = upper_M_VT,
                    interval_range = interval_range))

VT_res_G<-test_VT$VTadj[1:37] - UK_sim_G_avg$Virus_Tests
res_sum_sq_G_VT<-sum(VT_res_G^2)
sigma_VT_G<-sqrt((1/(37-5-1))*res_sum_sq_G_VT)
value_G_VT<-1.96*sigma_VT_G*sqrt(1+(1/37)+(UK_sim_G_avg$Virus_Tests-x_bar_VT)^2/((37-1)*sd_VT^2))


lower_G_VT<-UK_sim_G_avg$Virus_Tests-value_G_VT
upper_G_VT<-UK_sim_G_avg$Virus_Tests+value_G_VT

mean(interval_score(true_values = test_VT$VTadj[1:37],
                    lower = lower_G_VT,
                    upper = upper_G_VT,
                    interval_range = interval_range))

VT_res_I<-test_VT$VTadj[1:37] - UK_sim_I_avg$Virus_Tests
res_sum_sq_I_VT<-sum(VT_res_I^2)
sigma_VT_I<-sqrt((1/(37-5-1))*res_sum_sq_I_VT)
value_I_VT<-1.96*sigma_VT_I*sqrt(1+(1/37)+(UK_sim_I_avg$Virus_Tests-x_bar_VT)^2/((37-1)*sd_VT^2))

lower_I_VT<-UK_sim_I_avg$Virus_Tests-value_I_VT
upper_I_VT<-UK_sim_I_avg$Virus_Tests+value_I_VT

mean(interval_score(true_values = test_VT$VTadj[1:37],
                    lower = lower_I_VT,
                    upper = upper_I_VT,
                    interval_range = interval_range))
