# Çiftçi HÇ, Gümüş MG, Gümüş K. (2023). Hybrid VAR-VECM Model

library(tseries)
library(tidyverse)
library(stargazer)
library(readxl)
library(forecast)
library(seastests)
library(vars)
library(tsDyn)
library(urca)
library(BVAR)
library(RJDemetra)
library(Metrics)


# Reading meteorological data from Excel.

mois <- read_excel("C:/Users/Çağatay/Desktop/VAR Model/Nem.xlsx", col_names = TRUE)
temp <- read_excel("C:/Users/Çağatay/Desktop/VAR Model/Sıcaklık.xlsx", col_names = TRUE)

# Conversion of meteorological data into time series.

# Training sets
mois_train <- ts(mois, start = c(1950,1), end = c(2015,12), frequency = 12)
temp_train <- ts(temp, start = c(1950,1), end = c(2015,12), frequency = 12)

# Test sets
mois_test <- ts(mois, start = c(2016,1), end = c(2020,12), frequency = 12)
temp_test <- ts(temp, start = c(2016,1), end = c(2020,12), frequency = 12)

# TRAMO/SEATS with RJDemetra+

mois_tramo <- tramoseats(mois_train, spec = "RSAfull")
temp_tramo <- tramoseats(temp_train, spec = "RSAfull")

# Moisture components

mois_sa <- mois_tramo$decomposition$components[,2]
mois_t <- mois_tramo$decomposition$components[,3]
mois_s <- mois_tramo$decomposition$components[,4]
mois_i <- mois_tramo$decomposition$components[,5]

# Temperature components

temp_sa <- temp_tramo$decomposition$components[,2]
temp_t <- temp_tramo$decomposition$components[,3]
temp_s <- temp_tramo$decomposition$components[,4]
temp_i <- temp_tramo$decomposition$components[,5]

# Combining columns of meteorological time series.

y_sa <- cbind(mois_sa, temp_sa)
y_t <- cbind(mois_t,temp_t)
y_s <- cbind(mois_s,temp_s)
y_i <- cbind(mois_i,temp_i)

# Determining the appropriate lags according to various criteria.

lag_sa <- VARselect(y_sa, lag.max = 12, type = c("both"))
lag_y_t <- VARselect(y_t, lag.max = 12, type = c("both"))
lag_y_s <- VARselect(y_s, lag.max = 12, type = c("both"))
lag_y_i <- VARselect(y_i, lag.max = 12, type = c("both"))

################################################################
                        #VAR MODEL#
################################################################

# Estimation with appropriate VAR(p) lag.

estimate_y_sa <- VAR(y_sa, p=3, type = "both")
estimate_y_t <- VAR(y_t, p=12, type = "both")
estimate_y_s <- VAR(y_s, p=12, type = "both")
estimate_y_i <- VAR(y_i, p=8, type = "both")

# Forecast with VAR 

VAR_forecast_y_sa <- predict(estimate_y_sa, n.ahead = 60, ci = 0.95)
VAR_forecast_y_t <- predict(estimate_y_t, n.ahead = 60, ci = 0.95)
VAR_forecast_y_s <- predict(estimate_y_s, n.ahead = 60, ci = 0.95)
VAR_forecast_y_i <- predict(estimate_y_i, n.ahead = 60, ci = 0.95)


################################################################
                            #VECM MODEL#
################################################################

# Johansen Procedure for VAR

cajo_y_sa <- ca.jo(y_sa,type = c("eigen"), ecdet = c("trend"), K = 3,
                   spec=c("longrun"))
cajo_y_t <- ca.jo(y_t,type = c("eigen"), ecdet = c("trend"), K = 12,
                   spec=c("longrun"))
cajo_y_s <- ca.jo(y_s,type = c("eigen"), ecdet = c("trend"), K = 12,
                   spec=c("longrun"))
cajo_y_i <- ca.jo(y_i,type = c("eigen"), ecdet = c("trend"), K = 8,
                   spec=c("longrun"))


# VECM model building

VECM_y_sa <- VECM(y_sa, 3, r = 1,include = c("both"),beta = NULL, 
                  estim = c("2OLS"), LRinclude = c("both"),exogen = NULL)
VECM_y_t <- VECM(y_t, 12, r = 1,include = c("both"),beta = NULL, 
                  estim = c("2OLS"), LRinclude = c("both"),exogen = NULL)
VECM_y_s <- VECM(y_s, 12, r = 1,include = c("both"),beta = NULL, 
                  estim = c("2OLS"), LRinclude = c("both"),exogen = NULL)
VECM_y_i <- VECM(y_i, 8, r = 1,include = c("both"),beta = NULL, 
                  estim = c("2OLS"), LRinclude = c("both"),exogen = NULL)


# VECM to VAR

VECM_model_y_sa <- vec2var(cajo_y_sa, r = 1)
VECM_model_y_t <- vec2var(cajo_y_t, r = 1)
VECM_model_y_s <- vec2var(cajo_y_s, r = 1)
VECM_model_y_i <- vec2var(cajo_y_i, r = 1)

# Forecast with VECM

VECM_forecast_y_sa <- predict(VECM_model_y_sa, n.ahead = 120, ci = 0.95)
VECM_forecast_y_t <- predict(VECM_model_y_t, n.ahead = 120, ci = 0.95)
VECM_forecast_y_s <- predict(VECM_model_y_s, n.ahead = 120, ci = 0.95)
VECM_forecast_y_i <- predict(VECM_model_y_i, n.ahead = 120, ci = 0.95)

################################################################
                        #Hybrid VAR-VECM MODEL#
################################################################

hybrid_mois_sa <- (VAR_forecast_y_sa$fcst$mois_sa[,1]+VECM_forecast_y_sa$fcst$mois_sa[,1])/2
hybrid_mois_t <- (VAR_forecast_y_t$fcst$mois_t[,1]+VECM_forecast_y_t$fcst$mois_t[,1])/2
hybrid_mois_s <- (VAR_forecast_y_s$fcst$mois_s[,1]+VECM_forecast_y_s$fcst$mois_s[,1])/2
hybrid_mois_i <- (VAR_forecast_y_i$fcst$mois_i[,1]+VECM_forecast_y_i$fcst$mois_i[,1])/2


hybrid_temp_sa <-(VAR_forecast_y_sa$fcst$temp_sa[,1]+VECM_forecast_y_sa$fcst$temp_sa[,1])/2
hybrid_temp_t <- (VAR_forecast_y_t$fcst$temp_t[,1]+VECM_forecast_y_t$fcst$temp_t[,1])/2
hybrid_temp_s <- (VAR_forecast_y_s$fcst$temp_s[,1]+VECM_forecast_y_s$fcst$temp_s[,1])/2
hybrid_temp_i <- (VAR_forecast_y_i$fcst$temp_i[,1]+VECM_forecast_y_i$fcst$temp_i[,1])/2


# Hybrid VAR-VECM Forecast

hybrid_var_vecm_forecast_mois <- ((hybrid_mois_sa+hybrid_mois_t)/2)+hybrid_mois_s+hybrid_mois_i
hybrid_var_vecm_forecast_temp <- ((hybrid_temp_sa+hybrid_temp_t)/2)+hybrid_temp_s+hybrid_temp_i


hybrid_var_vecm_forecast_mois <- ts(hybrid_var_vecm_forecast_mois, 
                                    start = c(2016,1), end = c(2020,12), 
                                    frequency = 12)
hybrid_var_vecm_forecast_temp <- ts(hybrid_var_vecm_forecast_temp, 
                                    start = c(2016,1), end = c(2020,12), 
                                    frequency = 12)

RMSE_mois <- rmse(mois_test,hybrid_var_vecm_forecast_mois)
RMSE_temp <- rmse(temp_test,hybrid_var_vecm_forecast_temp)
