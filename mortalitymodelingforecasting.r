## Life table 

library(demography)
france<-hmd.mx("FRATNP","xxxxxx@hotmail.com","STK4500","FRATNP")

lt <- lifetable(france,type="period",years=1950)


lt

plot(lt)






## Fetching the data from the HMD.
## These are for the total population of France from 1816 to 2013.
## The two main quantities of interest are the Exposure (measured in person years) and the Death count.
## See "http://www.mortality.org/Public/ExplanatoryNotes.php#CompleteDataSeries" for more explanation.




country <- "FRATNP"
username <-"xxxxxx@hotmail.com"
password <- "STK4500"
path <- paste("http://www.mortality.org/hmd/", country, "/STATS/", "Exposures_1x1.txt", sep = "")
userpwd <- paste(username, ":", password, sep = "")
txt <- RCurl::getURL(path, userpwd = userpwd)
con <- textConnection(txt)
exposures <- try(read.table(con, skip = 2, header = TRUE, na.strings = "."),TRUE)

path <- paste("http://www.mortality.org/hmd/", country, "/STATS/", "Deaths_1x1.txt", sep = "")
userpwd <- paste(username, ":", password, sep = "")
txt <- RCurl::getURL(path, userpwd = userpwd)
con <- textConnection(txt)
deaths <- try(read.table(con, skip = 2, header = TRUE, na.strings = "."),TRUE)

EXPOSURE = exposures
DEATH = deaths

DEATH$Age=as.numeric(as.character(DEATH$Age))
DEATH$Age[is.na(DEATH$Age)]=110
EXPOSURE$Age=as.numeric(as.character(EXPOSURE$Age))
EXPOSURE$Age[is.na(EXPOSURE$Age)]=110


## Constructing the matrix of Mortality rates (compare with the mortality rates in HMD)

Mu<-DEATH[,3:5]/EXPOSURE[,3:5]
head(Mu)


## The force of mortality mu_{x,t} as a function of year and age

AGE<-as.numeric(DEATH$Age[1:111])
YEAR<-unique(as.numeric(DEATH$Year))
Mut<-matrix(Mu[,2], length(AGE), length(YEAR))
persp(AGE[1:90],YEAR,log(Mut[1:90,]),theta=-30,shade=T,col="light green",
      xlab="Age",ylab = "Year", zlab ="log(mu)",
      sub="Force of mortality as a function of time and age", d=3.5)

require(rgl)
persp3d(AGE[1:101],YEAR,log(Mut[1:101,]),theta=-30,shade=T,col="light green")



## The force of mortality mu_{x,t} as a function of age

library(rainbow)
rownames(Mut) <- AGE
colnames(Mut) <- YEAR
functionalTS <- fts(x = AGE[1:90], y = log(Mut[1:90,]), xname = "Age", yname = "Log Mortality rates")
plot(functionalTS, plot.type = "functions", plotlegend = TRUE, legendpos ="bottomright")


rownames(Mut) <- AGE
colnames(Mut) <- YEAR


MuM <- fts(x = AGE[1:90], y = log(Mut[1:90,80:198]), xname = "Ages", yname = "Log Mortality rates")


fboxplot(data = MuM, plot.type = "bivariate", type = "bag")
fboxplot(data = MuM, plot.type = "functional", type = "bag")

## Comparing $\mu_{x,t}$ for t=1900 and t=1990


year=1900

D1900=DEATH[DEATH$Year==year,]
E1900=EXPOSURE[EXPOSURE$Year==year,]
Mu1900 = D1900[,3:5]/E1900[,3:5]


year=1990

D1990=DEATH[DEATH$Year==year,]
E1990=EXPOSURE[EXPOSURE$Year==year,]
Mu1990 = D1990[,3:5]/E1990[,3:5]



par(mfrow=c(1,2))
plot(0:110,log(Mu1900[,1]), type ="o", col= "green", ylim = c(-10,1), xlab = "Age", ylab =  "Log Mortality rates",
     sub="Force of mortality 1900 as a function of age")
lines(log(Mu1900[,2]), type ="o", col= "red")

plot(0:110,log(Mu1990[,1]), type ="o", col= "green", ylim = c(-10,1), xlab = "Age", ylab =  "Log Mortality rates",
     sub="Force of mortality 1990 as a function of age")
lines(log(Mu1990[,2]), type ="o", col= "red")


## Comparing $p_{t}(x,x+1)$ for t=1900 and t=1990

PM1900=PF1900=matrix(0,111,111)
PM1990=PF1990=matrix(0,111,111)


for(x in 0:110){
  
  PM1900[x+1,1:(111-x)]=exp(-cumsum(Mu1900[(x+1):111,2]))
  PF1900[x+1,1:(111-x)]=exp(-cumsum(Mu1900[(x+1):111,1]))
  
  PM1990[x+1,1:(111-x)]=exp(-cumsum(Mu1990[(x+1):111,2]))
  PF1990[x+1,1:(111-x)]=exp(-cumsum(Mu1990[(x+1):111,1]))
}


par(mfrow=c(1,2))

plot(1:111,PF1900[1,],ylim=c(0,1),type="l",col="green", xlab = "Age", ylab =  "Survival probabilities at birth",
     sub=" 1900")
lines(PM1900[1,],type="l",col="red")

plot(1:111,PF1990[1,],ylim=c(0,1),type="l",col="green", xlab = "Age", ylab =  "Survival probabilities at birth",
     sub=" 1990")
lines(PM1990[1,],type="l",col="red")


## Estimation of the Lee-Carter model using $\textbf{demography}$ package


library(demography)
usa<-hmd.mx("USA","xxxxxx@hotmail.com","STK4500","USA")

usa.90<-extract.ages(usa,ages = 0:90 )

lc.male<- lca(usa.90, series = "male")
lc.female<- lca(usa.90, series = "female")

plot(lc.male)
plot(lc.female)


## Forcasting in the Lee-Carter model


forecast.lc.male <- forecast(lc.male, h=20)
forecast.lc.female <- forecast(lc.female, h=20)

plot(forecast.lc.male)
plot(forecast.lc.female)

kappat <- lc.female$kt
(fittedkt <- auto.arima(kappat))

plot(forecast(fittedkt, h=100))

## The LifeMetrics functions


### Fetching data from HMD

country <- "USA"
username <-"xxxxxx@hotmail.com"
password <- "STK4500"
path <- paste("http://www.mortality.org/hmd/", country, "/STATS/", "Exposures_1x1.txt", sep = "")
userpwd <- paste(username, ":", password, sep = "")
txt <- RCurl::getURL(path, userpwd = userpwd)
con <- textConnection(txt)
exposures <- try(read.table(con, skip = 2, header = TRUE, na.strings = "."),TRUE)

path <- paste("http://www.mortality.org/hmd/", country, "/STATS/", "Deaths_1x1.txt", sep = "")
userpwd <- paste(username, ":", password, sep = "")
txt <- RCurl::getURL(path, userpwd = userpwd)
con <- textConnection(txt)
deaths <- try(read.table(con, skip = 2, header = TRUE, na.strings = "."),TRUE)

deaths$Age=as.numeric(as.character(deaths$Age))
deaths$Age[is.na(deaths$Age)]=110
exposures$Age=as.numeric(as.character(exposures$Age))
exposures$Age[is.na(exposures$Age)]=110

### We extract ages up to 90 (Due to limitation in LifeMetrics)
DEATH_90<- deaths[deaths$Age<90,]
EXPOSURE_90<- exposures[exposures$Age<90,]

Ages <- unique(DEATH_90$Age)
Years <- unique(DEATH_90$Year)

ages <- seq(0,length(Ages)-1)
years <- seq(Years[1],Years[length(Years)])

n <- length(ages)
m <- length(years)
EXPOSUREM_tx <- t(matrix(EXPOSURE_90[,4],n,m))
EXPOSUREF_tx <- t(matrix(EXPOSURE_90[,3],n,m))
DEATHM_tx <- t(matrix(DEATH_90[,4],n,m))
DEATHF_tx <- t(matrix(DEATH_90[,3],n,m))

WeightFun <- matrix(1,m,n)


### Fitting the L-C model 

source("C:/Users/Amin/Documents/R/fitModels.r")


LeeCarterM <- fit701(xv=ages,yv=years,
                     etx=EXPOSUREM_tx,dtx=DEATHM_tx,
                     wa=WeightFun)
LeeCarterF <- fit701(xv=ages,yv=years,
                     etx=EXPOSUREF_tx,dtx=DEATHF_tx,
                     wa=WeightFun)


### Ploting the results

par(mfrow=c(2,2))

plot(LeeCarterM$x,LeeCarterM$beta1,type="l",col="red")
lines(LeeCarterF$x,LeeCarterF$beta1,col="green")

plot(LeeCarterM$x,LeeCarterM$beta2,type="l",col="red")
lines(LeeCarterF$x,LeeCarterF$beta2,col="green")

plot(LeeCarterM$y,LeeCarterM$kappa2,type="l",col="red", ylim = c(-100,60))
lines(LeeCarterF$y,LeeCarterF$kappa2,col="green")


### Forecasting with LifeMetrics


source("C:/Users/Amin/Documents/R/simModels.r")


forcastsLC <- sim2001(LeeCarterM$x,LeeCarterM$y,
                      LeeCarterM$beta1,LeeCarterM$beta2,
                      LeeCarterM$kappa2,
                      nsim = 10000,
                      nyears = 60,
                      tmax = 40,
                      x0 = 40)


par(mfrow=c(1,1))

fan(forcastsLC$xv,forcastsLC$tpa, pl = 1)



### Sensitivity of PV of an endowment with respect to longevity 

endowment <- array(0,length(forcastsLC$pa[1,]))

for(j in 1:length(endowment)){
  
  endowment[j] <- (1/(1+0.035)^40)*cumprod(forcastsLC$pa[,j])[length(forcastsLC$pa[,j])]
  
}

(sort(endowment)[length(endowment)*0.95])

sqrt(var(endowment))
mean(endowment)
hist(endowment,60)


par(mfrow=c(1,1))

plot(density(endowment))


### Effect of increased longevity on a life insurance annuity


LCM<-lca(usa, series = "male")




LCMf<-forecast(LCM,h=200)



A <- LCM$ax
B <- LCM$bx
K <- c(LCM$kt,LCM$kt[length(LCM$kt)]+LCMf$kt.f$mean) # Correct the length of k_t
Mu <- matrix(0,length(A),length(K))


for(i in 1:length(A)){
  
  for(j in 1:length(K)){
    
    Mu[i,j] <- exp(A[i]+B[i]*K[j])
    
  }
}

AGE <- seq(0,length(B)-1)
YEARS <- seq(years[1],years[1]+length(K)-1)
#persp(AGE, YEARS, log(Mu), theta=-30,shade=T,col="light green")


t <- 1960
x <- 40
s <- seq(0, 100-x-1)
MuD <- Mu[x+1+s, t+s-years[1]+2]
P_xt <- cumprod(exp(-diag(MuD)))


x <- 40
intRate <- 0.035
deferral <- 30
m <-x + deferral

horizon <- 120
valAnnuity <- rep(0,horizon + years[length(years)]-years[1]+1)

for(currYear in years[1]:(horizon + years[length(years)])){
  cohort <- seq(0,100-x-1)
  MuDD <- Mu[x+1+cohort,(cohort+currYear-years[1]+2)]
  P_xt <- cumprod(exp(-diag(MuDD)))
  Durations <- seq(0,deferral)
  val <- ((1/(1+intRate))^(m-x+Durations))*P_xt[m-x+Durations]
  valAnnuity[currYear-years[1]+1] <- sum(val)
  
}


plot(years[1]:(years[length(years)]+horizon),valAnnuity, xlab = "Year", ylab = "PV of an annuity")

