library(readr)
library(forecast)
library(tseries)
library(mgcv)

# 1.2
# reading data 
setwd("C:/Users/ruair/Downloads/Forecasting_Time_series")
condo_sales = read.csv("quebec_real_estate.csv",header=T)
str(condo_sales)


# plotting
values = condo_sales[4]/1000
values = ts(values, start=c(2014,1), frequency=12)
ts.plot(values, main="Montreal Median Price", ylab="Sale Price (1000's)", type="l")


# holdout
values.holdout <- window(values, start=2020,end=2021)
values <- window(values, start=1967,end=2019)
ts.plot(cbind(values, values.holdout), main="Quebec condo market  ",
        ylab="Price (1000's)", type="l", col=c("red", "blue"), lty=c(1, 2))


#smoothing the data
#Create equally spaced time points for fitting trends
time.pts = c(1:length(values))
time.pts = c(time.pts - min(time.pts))/max(time.pts)

# 1define mav/smoothing methods, and fit
values.mafilter.fit = filter(values, filter = rep(1/4, 4), sides = 2)
ma.fit = ma(values, order=2, centre=TRUE)
values.fit.ma = ts(ma.fit, start=c(2014, 1), frequency=12)

ksmooth.fit = ksmooth(time.pts, values, kernel = "box", bandwidth = 0.2)
values.fit.ksmooth = ts(ksmooth.fit$y,start=c(2014, 1),frequency=12)

loess.fit = loess(as.matrix(values)~time.pts, data=values, span=0.2)
values.fit.loess = ts(predict(loess.fit), start=c(2014, 1), frequency=12)

gam.fit = gam(values~s(time.pts))
values.fit.gam = ts(fitted(gam.fit),start=c(2014, 1),frequency=12)

# plot fits against values
lines(values.fit.ksmooth,lwd=1, lty=4 ,col="purple")
lines(values.mafilter.fit, col="red")
lines(values.fit.loess, col="orange", lty=4)
lines(values.fit.gam, col="violet")
lines(values.fit.ma, col="cyan", lwd=3)
legend(x="topleft",c("kernel smoothing", "filter", "loess", "gam", "mav"),lty = c(4, 1), col=c("purple"))
ts.plot(values,ylab="sale price", main="Observed Values vs smoothing methods")

values

# Decomposing the data
# use multiplicative because periods appear to be getting bigger
values.ma.decomp = decompose(values.fit.ma, type=c("multiplicative"))
plot(values.ma.decomp)

values.decomp.stl <- stl(na.omit(values.fit.ma[,1]), s.window="periodic")
plot(values.decomp.stl)

values.deseason = seasadj(values.decomp.stl)
plot(values.deseason, main="Deseasonalized time series")

print(values.decomp.stl)

# is the Data stationary?
values.fit.ma <- ts(na.omit(values.fit.ma), frequency=12, start=c(2014, 1))
acf(values.fit.ma)
adf.test(values.fit.ma)

# Differencing
values_diff12 = diff(values.fit.ma, lag = 12)
tm <- cbind(values, values_diff12)
plot(tm)

values_diff12_1 = diff(values_diff12)
tm <- cbind(values, values_diff12, values_diff12_1)
plot(tm)

acf(na.omit(values_diff12_1), main="ACF for series after 1st difference")

adf.test(na.omit(values_diff12_1))

diff.values <- na.omit(values_diff12_1)
par(mfrow=c(2,1))
acf(diff.values, main="")
pacf(diff.values, main="")

#  SARIMA(p, d, q)(P, D, Q) order
values.fit.sarima1 <- arima(values.fit.ma, order=c(1, 1, 0), seasonal = list(order = c(1, 1, 0), period = 12))
values.fit.sarima1
tsdisplay(residuals(values.fit.sarima1), lag.max=45, main='SARIMA Model 1 Residuals')
#forecast accuracy
forecast.sarima1 <- forecast(values.fit.sarima1, h=12) #12months
accuracy(forecast.sarima1)

values.fit.sarima2 <- arima(values.fit.ma, order=c(3, 1, 0), seasonal = list(order = c(1, 1, 0), period = 12))
values.fit.sarima2
#residuals
tsdisplay(residuals(values.fit.sarima2), lag.max=45, main='SARIMA Model 2 Residuals')
#forecast accuracy
forecast.sarima2 <- forecast(values.fit.sarima2, h=12) #12months
accuracy(forecast.sarima2)

values.fit.sarima3 <- arima(values.fit.ma, order=c(6, 1, 0), seasonal = list(order = c(1, 1, 0), period = 12))
values.fit.sarima3
tsdisplay(residuals(values.fit.sarima3), lag.max=45, main='SARIMA Model 3 Residuals')
forecast.sarima3 <- forecast(values.fit.sarima3, h=12) #12months
accuracy(forecast.sarima3)

# Auto ARIMA
values.fit.sarima4 <- auto.arima(ma.fit, seasonal = TRUE)
values.fit.sarima4
tsdisplay(residuals(values.fit.sarima4), lag.max=45, main='SARIMA Model 4 Residuals')
          
# holt winters
values.fit.hw1 <- ets(values.fit.ma, model="MAM", damped=FALSE)
values.fit.hw1
plot(values.fit.hw1)
accuracy(values.fit.hw1)

values.fit.hw2 <- ets(values.fit.ma, model="MAM",
                      damped=TRUE)
values.fit.hw2
plot(values.fit.hw2)
accuracy(values.fit.hw2)

# Forecast 
values = condo_sales[4]/1000
values = ts(values, start=c(2014,1), frequency=12)
ts.plot(values, main="Montreal Median Price", ylab="Sale Price (1000's)", type="l")

values.holdout <- window(values, start=2020,end=2021)
values2021 <- ts(values.holdout, start = c(2020, 1), frequency=12)
plot.ts(values2021)

# 6.2
forecast.values.sarima4 <- forecast(values.fit.sarima4, h=12)
plot(forecast.values.sarima4)
lines(values2021, col="cyan", lwd=2, lty=3)
accuracy(forecast.values.sarima4, values2021)

Box.test(residuals(values.fit.sarima4), type="Ljung-Box")

par(mfrow=c(2,2))
acf(na.omit(values.fit.sarima4$residuals))
pacf(na.omit(values.fit.sarima4$residuals))
plot(na.omit(values.fit.sarima4$residuals))
qqnorm(na.omit(values.fit.sarima4$residuals))
qqline(na.omit(values.fit.sarima4$residuals), col="cyan")

# Holt Winters
plot(forecast(values.fit.hw2))
lines(values2021, col="cyan", lwd=2, lty=3)
accuracy(forecast(values.fit.hw2), values2021)

Box.test(residuals(values.fit.hw2), type="Ljung-Box")

par(mfrow=c(2,2))
acf(values.fit.hw2$residuals)
pacf(values.fit.hw2$residuals)
plot(values.fit.hw2$residuals)
qqnorm(values.fit.hw2$residuals)
qqline(values.fit.hw2$residuals, col="cyan")

forecast.values.sarima4