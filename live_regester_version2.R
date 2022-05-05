library(readr)
library(forecast)
library(tseries)

# check the R version
R.Version()


# plotting
#setwd("C:/Users/sinea/Dropbox/second_stats/cso")
setwd("C:/Users/ruair/Downloads/Forecasting_Time_series")
live_register = read.csv("live_register_ireland_annual.csv",header=T)
values = live_register[8]/1000
values = ts(values, start=1967, end=2021,frequency=1)
ts.plot(values, main="Live Register Monthly Mean", ylab="No. of People (1000's)", type="l")

#output the values 
values

#check for missing values
complete <- TRUE
for(c in complete.cases(values)) {
  if(!c){
    complete == FALSE
  }
}
if(complete){
  print("No missing values were found")
} else {
  print ("Missing values found")
}

#set up a holdout set of observations for comparison with the forecast
values.holdout <- window(values, start=2017,end=2021)
values <- window(values, start=1967,end=2016)

#plot the shortened series and the holdout
ts.plot(cbind(values, values.holdout), main="Live Register Monthly Mean",
        ylab="No. of People (1000's)", type="l", col=c("red", "blue"), lty=c(1, 2))


# Check if the data is stationary using the KPSS and ADF tests
adf.test(values)

kpss.test(values, null="Trend")
kpss.test(values, null="Level")


# transform, difference 
values_diff1 = diff(values, lag = 1)
values_diff2 = diff(values_diff1)
tdiff <- cbind(values, values_diff1, values_diff2)
plot(tdiff, main="1st and 2nd order Differencing with Lag=1")

#Get ADF and ACF for 1st Difference

par(mfrow=c(2,1)) 
acf(values_diff1) 
adf.test(values_diff1)

#Get ADF and ACF for 2nd Difference
par(mfrow=c(2,1)) 
acf(values_diff2) 
adf.test(values_diff2)


# Decomposition, fitting moving average

#Create equally spaced time points for fitting trends
time.pts = c(1:length(values))
time.pts = c(time.pts - min(time.pts))/max(time.pts)

#get the moving average 
ma.values <- ma(values, order=4, centre = FALSE)

#define mav method and fit
mav.fit = ksmooth(time.pts, values, kernel = "box", bandwidth = 1.1)
values.fit.mav = ts(mav.fit$y,start=1967,frequency=1)

# plot mav.fit against values
ts.plot(values,ylab="number of people", main="Observed Values vs MA vs Kernel smoothing")
lines(ma.values,lwd=2, lty=4, col="violet")
lines(values.fit.mav, col="red")
legend(x="topleft", c("Observed", "MA", "K smoothing"), col=c("gray10","violet", "red"), lty=c(1, 4))

# 3 detrending
loess_fit25 <- loess(as.matrix(values)~time.pts, data=values, span=0.25) #25%smoothing span
loess_fit50 <- loess(as.matrix(values)~time.pts, data=values, span=0.5) #50%smoothing span
smoothed.values.loess25 <- ts(fitted(loess_fit25), start=1967)
smoothed.values.loess50 <- ts(predict(loess_fit50), start=1967)
plot(cbind(ma.values, smoothed.values.loess50), main = "De-Trending")


# LOESS
loess.fit = loess(as.matrix(values)~time.pts, degree=2)
values.fit.loess = ts(fitted(loess.fit),start=1967)
plot(values.fit.loess)

# plot LOESS against observations and MA
ts.plot(values,ylab="Number of People", main="Observations, MA & LOESS")
lines(ma.values,lwd=2, lty=4, col="blue")
lines(values.fit.loess, col="red")
legend(x="topleft", c("Observed", "MA", "LOESS"), col=c("gray10","violet", "red"), lty=c(1, 4))


# plot LOESS against observ and ma
ts.plot(values,ylab="# of PAtent Applications", main="Observed Values vs MA vs LOESS")
lines(ma.values,lwd=2, lty=4, col="violet")
lines(smoothed.values.loess50, col="red")
legend(x="topleft", c("Observed", "MA (order 4)", "LOESS (span 0.5)"), col=c("gray10","violet", "red"))

values.rand.ma <- values-ma.values
values.rand.loess <- ts(loess_fit50$residuals, start = 1977)
plot(values.rand.loess, col="tomato", pty="l", main="Random component", ylim=c(1000, -1000))
lines(values.rand.ma, lty=2)
legend(x="topleft", legend=c("Values - MA estimated trend", "Values - LOESS estimated trend"),
       col=c("gray15", "tomato"), lty=c(2, 1))



# 4.1 ETS
HoltWinters(values_diff2, beta=FALSE, gamma=FALSE)
plot(HoltWinters(values_diff2, beta=FALSE, gamma=FALSE))


values.hsmooth <- HoltWinters(values, gamma=FALSE)
values.hsmooth
plot(values.hsmooth)


# 4.2 ARIMA
arima(values_diff2)

acf(values_diff2)
pacf(values_diff2)


# Model comparison
fets <- function(x, h) {
  forecast(ets(x), h = h)
}
farima <- function(x, h) {
  forecast(auto.arima(x), h=h)
}
# Compute CV errors for ETS
res.ets <- tsCV(values, fets, h=5)

# Compute CV errors for ARIMA
res.autoarima <- tsCV(values, fautoarima, h=5)

# Find MSE of each model class
mean(res.ets^2, na.rm=TRUE)

mean(res.autoarima^2, na.rm=TRUE)

fets <- function(x, h) {
  forecast(ets(x, model="MAN"), h = h)
}
farima <- function(x, h, order) {
  forecast(arima(x, order = c(0,2,2)), h = h)
}
# Compute CV errors for ETS
res.MANets <- tsCV(values, fets, h=1)

# Compute CV errors for ARIMA
res.arima0_2_2 <- tsCV(values, farima, h=1)

# Find MSE of each model class
sqrt(mean(res.MANets^2, na.rm=TRUE))

sqrt(mean(res.arima0_2_2^2, na.rm=TRUE))

res.MANets

res.arima0_2_2

# Forecasting
forecast.fit.MANets <- forecast(values.fit.MANets, h=3, bootstrap = TRUE)
plot(forecast.fit.MANets)
lines(values.holdout, col="tomato", lwd=2)

accuracy(forecast.fit.MANets, values.holdout)

forecast.fit.arima0_2_2 <- forecast(values.fit.arima0_2_2, h=3, bootstrap = TRUE)
plot(forecast.fit.arima0_2_2)
lines(values.holdout, col="tomato", lwd=2)

accuracy(forecast.fit.arima0_2_2, values.holdout)



