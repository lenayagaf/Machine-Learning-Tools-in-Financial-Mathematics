# Clean the global environment
rm(list = ls())

# Set the working directory
dir_data<-".MSc-Thesis/Return volatility estimation and forecast in R/"
setwd(dir_data)

# Load the required libraries
library(data.table)
library(estimatr)
library(ggplot2)
library(patchwork)
library(stringr)
library(FinTS)
library(xtable)
library(forecast)
library(tseries)
library(moments)
library(dplyr)
library(fGarch)
library(stats)

# Read the data
data_dummy<-fread("eur_gbp_data.csv")

# Number of observations in the training set (70% of the data)
n_train <- round(dim(data_dummy)[1]*0.7)

# Filter the training dataset
data_dummy<-data_dummy[1:n_train,]

# Standardize the data (mean = 1.637077e-06, standard deviation = 0.0002286759)
data_dummy$Return1 <- scale(data_dummy$Return, center = TRUE, scale = TRUE)

# Calculate summary statistics, skewness, kurtosis and Jarque-Bera test statistic for the non-standardized data
summary(data_dummy$Return)
sd(data_dummy$Return)
skewness(data_dummy$Return)
kurtosis(data_dummy$Return)
jarque.bera.test(data_dummy$Return)

# Calculate summary statistics, skewness, kurtosis and Jarque-Bera test statistic for the standardized data
summary(data_dummy$Return1)
sd(data_dummy$Return1)
jarque.bera.test(data_dummy$Return1)
skewness(data_dummy$Return1)

# Empirical and simulated normal density for Standardized EUR/GBP 5-minute Log Return
ggplot(data_dummy, aes(x=Return1))+geom_density(color="red", linewidth = 1) +
  theme(axis.text=element_text(size=50), axis.title=element_text(size=50))+
  labs(y = "", x = "") + stat_function(fun = dnorm, args = list(mean = mean(data_dummy$Return1), sd=sd(data_dummy$Return1)), color="green", linewidth = 1) +
  xlim(-10, 10) + theme_classic() + ggtitle("Density of Standardized EUR/GBP 5-minute Log Returns") +
  theme(plot.title = element_text(size = 15, face = "bold"))

# Plot a normal Quantile-Quantile plot for Standardized EUR/GBP 5-minute Log Return
qqnorm(data_dummy$Return1, main = "QQ plot for Standardized EUR/GBP 5-minute Log Returns", pch = 19, col = "blue", cex = 0.5,
       xlab = "Theoretical Quantiles",
       ylab = "Quantiles of Standardized EUR/GBP 5-minute Log Returns", datax=F)
qqline(data_dummy$Return1, lwd = 3, col="orange")

# Non-robust pure diurnal dummy variable model for Standardized EUR/GBP 5-minute Log Returns
model_return_non_robust<-lm(Return1~D1+D2+D3+D4+D5+D6+D7+D8+D9+D10+D11+D12+D13+D14+D15+D16-1, data=data_dummy)
residuals_model_return_non_robust <- residuals(model_return_non_robust)

    # Residual analysis of non-robust pure diurnal dummy variable model for Standardized EUR/GBP 5-minute Log Returns
    data_dummy1 <- cbind(data_dummy, residuals_model_return_non_robust)
    ggplot(data=data_dummy1, aes(x=Index, y=residuals_model_return_non_robust^2)) + geom_point() + theme_classic() +
      theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) + 
      labs(y = "", x = "") + ggtitle("Squared Seasonally Adjusted Log Return Series") +
      theme(plot.title = element_text(size = 20, face = "bold"))
    
    # Autocorrelation plot 
    Acf(residuals_model_return_non_robust, main="Non-robust Seasonally Adjusted Log Return Series ACF")
    # Arch-LM test 
    ArchTest(residuals_model_return_non_robust)

# Robust pure diurnal dummy variable model for Standardized EUR/GBP 5-minute Log Returns
model_return_robust<-lm_robust(Return1~D1+D2+D3+D4+D5+D6+D7+D8+D9+D10+D11+D12+D13+D14+D15+D16-1, data=data_dummy, se_type = "HC3")
#format(model_return_robust, scientific = FALSE)
residuals_model_return_robust <- data_dummy$Return1 - model_return_robust$fitted.values

Acf(residuals_model_return_robust)
kurtosis(residuals_model_return_robust)

residuals_model_return_robust<-as.data.frame(residuals_model_return_robust)
colnames(residuals_model_return_robust)<-c("Return")

robust_model_LB<-lm(return~lag(return, 1)+lag(return, 2)+
     lag(return, 3)+lag(return, 4)+lag(return, 5), data=residuals_model_return_robust)

mean_est<-0
stat_chi<-((robust_model_LB$coefficients[2]-mean_est)/(sqrt(diag(vcov(robust_model_LB)))[2]))^2
for (i in 3:6) {
  stat_chi<-stat_chi+((robust_model_LB$coefficients[i]-mean_est)/sqrt(diag(vcov(robust_model_LB)))[i])^2
}
stat_chi

# Plot a specific interval for visual clarity
# Choose length of the interval to be plotted
interval_length<-500
Standardized_Return <- rep("Standardized Return", interval_length)
Deseasonalized_Standardized_Return <- rep("Deseasonalized & Standardized Return", interval_length)
# Choose a starting point
l1 <- 18001
l2 <- l1+interval_length - 1

df_returns_deseason1 <- cbind(data_dummy$Index[l1:l2], data_dummy$Return1[l1:l2], Standardized_Return)
df_returns_deseason2 <- cbind(data_dummy$Index[l1:l2], residuals_model_return_robust$Return[l1:l2], Deseasonalized_Standardized_Return)

df_returns_deseason<-as.data.frame(rbind(df_returns_deseason1, df_returns_deseason2))
colnames(df_returns_deseason)<-c("Index", "Return", "Type")
df_returns_deseason$Index<-as.numeric(df_returns_deseason$Index)
df_returns_deseason$Return<-as.numeric(df_returns_deseason$Return)

ggplot(df_returns_deseason, aes(x=Index, group = factor(Type), colour=factor(Type), y=Return)) + geom_line(linewidth=1) + 
  theme(axis.text=element_text(size=50), axis.title=element_text(size=50)) +
  labs(y = "", x = "") +
  theme_classic() + ggtitle("Standardized and Deseasonalized & Standardized Return Series") +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  scale_color_manual(name = "", breaks = c("Standardized Return", "Deseasonalized & Standardized Return"), values = c("green", "blue")) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
  theme(legend.title = element_text(size=15), legend.text = element_text(size=15))

mean_est<-mean(model_return_robust$coefficients)
stat_chi<-((model_return_robust$coefficients[1]-mean_est)/(model_return_robust$std.error[1]))^2
for (i in 2:16) {
  stat_chi<-stat_chi+((model_return_robust$coefficients[i]-mean_est)/model_return_robust$std.error[i])^2
}

# Plot autocorrelation function 
Acf(residuals_model_return_robust, main="EUR/GBP 5-minute Deseasonalized Log Return ACF")
# Plot partial autocorrelation function
pacf(residuals_model_return_robust, main="EUR/GBP 5-minute Deseasonalized Log Returns PACF")

# ARMA model estimations 
ma1<-Arima(residuals_model_return_robust, order=c(0, 0, 1), include.mean = F)
ma2<-Arima(residuals_model_return_robust, order=c(0, 0, 2), include.mean = F)
ma3<-Arima(residuals_model_return_robust, order=c(0, 0, 3), include.mean = F)
ar1<-Arima(residuals_model_return_robust, order=c(1, 0, 0), include.mean = F)
ar2<-Arima(residuals_model_return_robust, order=c(2, 0, 0), include.mean = F)
ar3<-Arima(residuals_model_return_robust, order=c(3, 0, 0), include.mean = F)
arima1<-Arima(residuals_model_return_robust, order=c(1, 0, 1), include.mean = F)
arima2<-Arima(residuals_model_return_robust, order=c(1, 0, 2), include.mean = F)
arima3<-Arima(residuals_model_return_robust, order=c(1, 0, 3), include.mean = F)
  arima12<-Arima(residuals_model_return_robust, order=c(2, 0, 1), include.mean = F)
arima22<-Arima(residuals_model_return_robust, order=c(2, 0, 2), include.mean = F)
arima32<-Arima(residuals_model_return_robust, order=c(2, 0, 3), include.mean = F)
arima123<-Arima(residuals_model_return_robust, order=c(3, 0, 1), include.mean = F)
arima223<-Arima(residuals_model_return_robust, order=c(3, 0, 2), include.mean = F)
arima323<-Arima(residuals_model_return_robust, order=c(3, 0, 3), include.mean = F)
arima41<-Arima(residuals_model_return_robust, order=c(4, 0, 1), include.mean = F)
arima42<-Arima(residuals_model_return_robust, order=c(4, 0, 2), include.mean = F)
arima43<-Arima(residuals_model_return_robust, order=c(4, 0, 3), include.mean = F)

summary(arima12)

Acf(arima43$residuals)
Acf(arima12$residuals)

vec_aic <- c(ma1$aic, ma2$aic, ma3$aic, ar1$aic, ar2$aic, ar3$aic, arima1$aic, arima2$aic, arima3$aic, arima12$aic, arima22$aic, arima32$aic, arima123$aic, arima223$aic, arima323$aic, arima41$aic, arima42$aic, arima43$aic)
vec_aicc <- c(arima1$aicc, arima2$aicc, arima3$aicc, arima12$aicc, arima22$aicc, arima32$aicc, arima123$aicc, arima223$aicc, arima323$aicc)
vec_bic <- c(ma1$bic, ma2$bic, ma3$bic, ar1$bic, ar2$bic, ar3$bic, arima1$bic, arima2$bic, arima3$bic, arima12$bic, arima22$bic, arima32$bic, arima123$bic, arima223$bic, arima323$bic, arima41$bic, arima42$bic, arima43$bic)
model<-c("ARMA(0,1)", "ARMA(0,2)", "ARMA(0,3)", "ARMA(1,0)", "ARMA(2,0)", "ARMA(3,0)", "ARMA(1, 1)", "ARMA(1, 2)", "ARMA(1, 3)", "ARMA(2, 1)",  "ARMA(2, 2)", "ARMA(2, 3)", "ARMA(3, 1)", "ARMA(3, 2)", "ARMA(3, 3)", "ARMA(4, 1)", "ARMA(4, 2)", "ARMA(4, 3)")
aic <- rep("AIC", 18)
bic <- rep("BIC", 18)

cbind(aic, vec_aic)
cbind(bic, vec_bic)

df_ic <- as.data.frame(rbind(cbind(aic, vec_aic, model),cbind(bic, vec_bic, model)))
colnames(df_ic)<-c("ic", "value", "model")
df_ic$value <- as.numeric(df_ic$value)

df_ic_aic<-filter(df_ic, ic == "AIC")
df_ic_bic<-filter(df_ic, ic == "BIC")

minaic <- df_ic_aic[which(df_ic_aic$value == min(df_ic_aic$value)),]
minbic <- df_ic_bic[which(df_ic_bic$value == min(df_ic_bic$value)),]

ggplot(df_ic, aes(x=factor(model), group = factor(ic), colour=factor(ic), y=value)) + geom_line(linewidth=1) + 
  geom_point(data=minaic, aes(x=factor(model), y = value), color = "red") +
  geom_point(data=minbic, aes(x=factor(model), y = value), color = "red") +
  theme(axis.text=element_text(size=50), axis.title=element_text(size=50)) +
  labs(y = "", x = "") +
  theme_classic() + ggtitle("Information Criteria by ARMA(p,q) Models") +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  scale_color_manual(name = "Information criterion", breaks = c("AIC", "BIC"), values = c("green", "blue")) +
  scale_x_discrete("", labels = factor(model), breaks = factor(model)) + theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
  theme(legend.title = element_text(size=15), legend.text = element_text(size=15))

# Automatic SARMA model selection (self-check)
auto.arima(residuals_model_return_robust,max.p = 5,max.q = 5,max.P = 5,max.Q = 5,max.d = 3,seasonal = TRUE,ic = 'aicc')

data_dummy$Residuals <- arima12$residuals

data_dummy$Residuals <- as.vector(data_dummy$Residuals)

ggplot(data=data_dummy, aes(x=Index, y=Residuals^2)) + geom_point() + theme_classic() +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) + 
  labs(y = "", x = "") + ggtitle("Fitted Residuals of ARMA(2,1) Model") + theme(plot.title = element_text(size = 15, face = "bold"))

# ARMA Residual diagnostics
ArchTest(arima12$residuals)

Acf(arima12$residuals, main="Autocorrelation Function for Fitted Residuals")
pacf(arima12$residuals, main="Autocorrelation Function for Fitted Residuals")

model_residuals <- lm_robust(arima12$residuals~1, se_type = "HC3")
format(model_residuals, scientific = FALSE)

Box.test(arima12$residuals, type = "Ljung")
Box.test(arima43$residuals, type = "Ljung")

robust_model_LB10<-lm_robust(Residuals~lag(Residuals, 1)+lag(Residuals, 2)+
                       lag(Residuals, 3)+lag(Residuals, 4)+lag(Residuals, 5)+
                       lag(Residuals, 6)+lag(Residuals, 7)+
                       lag(Residuals, 8)+lag(Residuals, 9)+lag(Residuals, 10), data=data_dummy, se_type = "HC3")

mean_est<-0
stat_chi<-((robust_model_LB10$coefficients[2]-mean_est)/robust_model_LB10$std.error[2])^2
for (i in 3:11) {
  stat_chi<-stat_chi+((robust_model_LB10$coefficients[i]-mean_est)/robust_model_LB10$std.error[i])^2
}
Box.test(data_dummy$Residuals, lag = 10, type = "Ljung")

robust_model_LB20<-lm_robust(Residuals~lag(Residuals, 1)+lag(Residuals, 2)+
                       lag(Residuals, 3)+lag(Residuals, 4)+lag(Residuals, 5)+
                       lag(Residuals, 6)+lag(Residuals, 7)+
                       lag(Residuals, 8)+lag(Residuals, 9)+lag(Residuals, 10)+
                       lag(Residuals, 11)+lag(Residuals, 12)+
                       lag(Residuals, 13)+lag(Residuals, 14)+lag(Residuals, 15)+
                       lag(Residuals, 16)+lag(Residuals, 17)+
                       lag(Residuals, 18)+lag(Residuals, 19)+lag(Residuals, 20), data=data_dummy, se_type = "HC3")

mean_est<-0
stat_chi<-((robust_model_LB20$coefficients[2]-mean_est)/robust_model_LB20$std.error[2])^2
for (i in 3:21) {
  stat_chi<-stat_chi+((robust_model_LB20$coefficients[i]-mean_est)/robust_model_LB20$std.error[i])^2
}
Box.test(data_dummy$Residuals, lag = 20, type = "Ljung")

robust_model_LB40<-lm_robust(Residuals~lag(Residuals, 1)+lag(Residuals, 2)+
                               lag(Residuals, 3)+lag(Residuals, 4)+lag(Residuals, 5)+
                               lag(Residuals, 6)+lag(Residuals, 7)+
                               lag(Residuals, 8)+lag(Residuals, 9)+lag(Residuals, 10)+
                               lag(Residuals, 11)+lag(Residuals, 12)+
                               lag(Residuals, 13)+lag(Residuals, 14)+lag(Residuals, 15)+
                               lag(Residuals, 16)+lag(Residuals, 17)+
                               lag(Residuals, 18)+lag(Residuals, 19)+lag(Residuals, 20)+
                               lag(Residuals, 21)+lag(Residuals, 22)+
                               lag(Residuals, 23)+lag(Residuals, 24)+lag(Residuals, 25)+
                               lag(Residuals, 26)+lag(Residuals, 27)+
                               lag(Residuals, 28)+lag(Residuals, 29)+lag(Residuals, 30)+
                               lag(Residuals, 31)+lag(Residuals, 32)+
                               lag(Residuals, 33)+lag(Residuals, 34)+lag(Residuals, 35)+
                               lag(Residuals, 36)+lag(Residuals, 37)+
                               lag(Residuals, 38)+lag(Residuals, 39)+lag(Residuals, 40), data=data_dummy, se_type = "HC3")

mean_est<-0
stat_chi<-((robust_model_LB40$coefficients[2]-mean_est)/robust_model_LB40$std.error[2])^2
for (i in 3:41) {
  stat_chi<-stat_chi+((robust_model_LB40$coefficients[i]-mean_est)/robust_model_LB40$std.error[i])^2
}
Box.test(data_dummy$Residuals, lag = 40, type = "Ljung")


T <- length(arima12$residuals)
k <- 3
se <- sum((arima12$residuals - mean(arima12$residuals))^2)/(T-k-1)

mean(arima12$residuals) * sqrt(T) / se

#Volatility 

arima12$residuals

data_dummy$Residuals <- arima12$residuals
log(data_dummy$Residuals^2)

model_volatility_robust<-lm_robust(log(Residuals^2)~D1+D2+D3+D4+D5+D6+D7+D8+D9+D10+D11+D12+D13+D14+D15+D16-1, data=data_dummy, se_type = "HC3")

mean_est<-mean(model_volatility_robust$coefficients)
stat_chi<-((model_volatility_robust$coefficients[1]-mean_est)/(model_volatility_robust$std.error[1]))^2
for (i in 2:16) {
  stat_chi<-stat_chi+((model_volatility_robust$coefficients[i]-mean_est)/model_volatility_robust$std.error[i])^2
}

model_volatility_robust_residuals <- log(data_dummy$Residuals^2) - model_volatility_robust$fitted.values
  
S_squared_epsilon <- exp(model_volatility_robust$fitted.values) * sum(exp(model_volatility_robust_residuals)) / length(model_volatility_robust$fitted.values)

epsilon_star <- data_dummy$Residuals / sqrt(S_squared_epsilon)

mean(epsilon_star)
sd(epsilon_star)
skewness(epsilon_star)
kurtosis(epsilon_star)

num<-seq(1, length(epsilon_star), 1)
df <- as.data.frame(cbind(num, epsilon_star))

# epsilon_star

ggplot(df, aes(x=epsilon_star))+geom_density(color="red", linewidth = 1) +
  theme(axis.text=element_text(size=50), axis.title=element_text(size=50))+
  labs(y = "", x = "") + stat_function(fun = dged, args = list(mean = mean(df$epsilon_star), sd=sd(df$epsilon_star)), color="green", linewidth = 1) + 
  xlim(-10, 10) + theme_classic() + ggtitle("Density of Standardized EUR/GBP 5-minute Log Returns") +
  theme(plot.title = element_text(size = 15, face = "bold"))

data_dummy$epsilon_star <- epsilon_star

#QQ plot
N         <- length(epsilon_star);
PERCS     <- ((1:N)-0.5)/N;
QUANTILES1 <- qsged(PERCS, mean = mean(epsilon_star), sd = sd(epsilon_star), xi = skewness(epsilon_star), nu = 3);
QUANTILES2 <- qsged(PERCS, mean = mean(epsilon_star), sd = sd(epsilon_star), xi = skewness(epsilon_star), nu = 5);
QUANTILES3 <- qged(PERCS, mean = mean(epsilon_star), sd = sd(epsilon_star), nu = 3)
QUANTILES4 <- qged(PERCS, mean = mean(epsilon_star), sd = sd(epsilon_star), nu = 5)
QUANTILES5 <- qsstd(PERCS, mean = mean(epsilon_star), sd = sd(epsilon_star), xi = skewness(epsilon_star), nu = 3)
QUANTILES6 <- qsstd(PERCS, mean = mean(epsilon_star), sd = sd(epsilon_star), xi = skewness(epsilon_star), nu = 5)
QUANTILES7 <- qstd(PERCS, mean = mean(epsilon_star), sd = sd(epsilon_star), nu = 3)
QUANTILES8 <- qstd(PERCS, mean = mean(epsilon_star), sd = sd(epsilon_star), nu = 5)


PLOTDATA1 <- data.frame(Sample = sort(epsilon_star),
                       Theoretical = QUANTILES1)
PLOTDATA2 <- data.frame(Sample = sort(epsilon_star),
                       Theoretical = QUANTILES2)
PLOTDATA3 <- data.frame(Sample = sort(epsilon_star),
                       Theoretical = QUANTILES3)
PLOTDATA4 <- data.frame(Sample = sort(epsilon_star),
                       Theoretical = QUANTILES4)
PLOTDATA5 <- data.frame(Sample = sort(epsilon_star),
                        Theoretical = QUANTILES5)
PLOTDATA6 <- data.frame(Sample = sort(epsilon_star),
                        Theoretical = QUANTILES6)
PLOTDATA7 <- data.frame(Sample = sort(epsilon_star),
                        Theoretical = QUANTILES7)
PLOTDATA8 <- data.frame(Sample = sort(epsilon_star),
                        Theoretical = QUANTILES8)

#Generate custom QQ plot
theme_update(plot.title    = element_text(size = 25, hjust = 0.5),
             plot.subtitle = element_text(size = 1, hjust = 0.5),
             axis.title.x  = element_text(size = 10, hjust = 0.5),
             axis.title.y  = element_text(size = 10, vjust = 0.5));

p1 <- ggplot(data = PLOTDATA1, aes(x = Theoretical, y = Sample)) +
  geom_point(size = 1, colour = 'blue') + 
  geom_abline(intercept = 0, slope = 1, col = "orange", linewidth = 1) + 
  xlab('Quantiles of Skew generalized error distribution (df=3)') + 
  ylab('Quantiles of Standardized GARCH shock') + theme_classic()

p2 <- ggplot(data = PLOTDATA2, aes(x = Theoretical, y = Sample)) +
  geom_point(size = 1, colour = 'blue') + 
  geom_abline(intercept = 0, slope = 1, col = "orange", linewidth = 1) + 
  xlab('Quantiles of Skew generalized error distribution (df=5)') + 
  ylab('Quantiles of Standardized GARCH shock') + theme_classic()

p3 <- ggplot(data = PLOTDATA3, aes(x = Theoretical, y = Sample)) +
  geom_point(size = 1, colour = 'blue') + 
  geom_abline(intercept = 0, slope = 1, col = "orange", linewidth = 1) + 
  xlab('Quantiles of Standardized generalized error distribution (df=3)') + 
  ylab('Quantiles of Standardized GARCH shock') + theme_classic()

p4 <- ggplot(data = PLOTDATA4, aes(x = Theoretical, y = Sample)) +
  geom_point(size = 1, colour = 'blue') + 
  geom_abline(intercept = 0, slope = 1, col = "orange", linewidth = 1) + 
  xlab('Quantiles of Standardized generalized error distribution (df=5)') + 
  ylab('Quantiles of Standardized GARCH shock') + theme_classic()

p5 <- ggplot(data = PLOTDATA5, aes(x = Theoretical, y = Sample)) +
  geom_point(size = 1, colour = 'blue') + 
  geom_abline(intercept = 0, slope = 1, col = "orange", linewidth = 1) + 
  xlab('Quantiles of  Skew Student-t distribution (df=3)') + 
  ylab('Quantiles of Standardized GARCH shock') + theme_classic()

p6 <- ggplot(data = PLOTDATA6, aes(x = Theoretical, y = Sample)) +
  geom_point(size = 1, colour = 'blue') + 
  geom_abline(intercept = 0, slope = 1, col = "orange", linewidth = 1) + 
  xlab('Quantiles of  Skew Student-t distribution (df=5)') + 
  ylab('Quantiles of Standardized GARCH shock') + theme_classic()

p7 <- ggplot(data = PLOTDATA7, aes(x = Theoretical, y = Sample)) +
  geom_point(size = 1, colour = 'blue') + 
  geom_abline(intercept = 0, slope = 1, col = "orange", linewidth = 1) + 
  xlab('Quantiles of Standardized Student-t distribution (df=3)') + 
  ylab('Quantiles of Standardized GARCH shock') + theme_classic()

p8 <- ggplot(data = PLOTDATA8, aes(x = Theoretical, y = Sample)) +
  geom_point(size = 1, colour = 'blue') + 
  geom_abline(intercept = 0, slope = 1, col = "orange", linewidth = 1) + 
  xlab('Quantiles of Standardized Student-t distribution (df=5)') + 
  ylab('Quantiles of Standardized GARCH shock') + theme_classic()

p1+p2
p3+p4
p5+p6
p7+p8

# Fit GARCH(1,1) model with a pre-specified conditional distribution for deseasonalized residual epsilon_star
m1 <- garchFit(formula = ~ garch(1, 1), epsilon_star, cond.dist = "std", include.shape = F, shape=3, trace=F)
coef(m1)["alpha1"] + coef(m1)["beta1"]
as.data.frame(summary(m1))
options(scipen=999)
xtable(format(m1@fit$matcoef, digits=1), type = "latex")
#AIC      BIC      SIC     HQIC 
#2.168601 2.169653 2.168601 2.168938 

m2 <- garchFit(formula = ~ garch(1, 1), epsilon_star, cond.dist = "std", include.shape = F, shape=5, trace=F)
coef(m2)["alpha1"] + coef(m2)["beta1"]
summary(m2)
#AIC      BIC      SIC     HQIC 
#2.142785 2.143837 2.142785 2.143122 
xtable(format(m2@fit$matcoef, digits=1), type = "latex")

deseasonalized_volatility <- volatility(m2)
true_volatility <- sqrt(S_squared_epsilon) * deseasonalized_volatility
index <- seq(1, length(true_volatility), 1)

df_volatility <- as.data.frame(cbind(index, deseasonalized_volatility, true_volatility))

# Conduct a volatility forecast optimality test via alternative Mincer-Zarnowitz regression 
mincer_zarnowitz <- as.data.frame(cbind(index, arima12$residuals^2/true_volatility^2, 1/true_volatility^2))
colnames(mincer_zarnowitz)<-c("i", "y", "x")
mincer_zarnowitz_model <- lm(y~x, data=mincer_zarnowitz)
xtable(mincer_zarnowitz_model)

# Define the interval in terms of index to plot
# Choose interval length
interval_length<-500
# Choose a starting point
l1 <-21000  
l2 <-l1+interval_length - 1

# Plot deseasonalized and true volatility for a pre-specified interval (Deseasonalized and True Volatility of Standardized EUR/GBP 5-minute Log Returns)
ggplot(df_volatility[l1:l2,], aes(x=index)) + 
  geom_line(aes(y=deseasonalized_volatility, col = "Deseasonalized volatility", linetype = "Deseasonalized volatility")) +
  geom_line(aes(y=true_volatility, col = "True volatility", linetype = "True volatility")) +
  scale_color_manual(name = "", breaks = c("Deseasonalized volatility", "True volatility"), values = c("steelblue", "darkred")) +
  scale_linetype_manual(name = "", breaks = c("Deseasonalized volatility", "True volatility"), values = c("dashed", "solid")) +
  labs(y = "", x = "") +
  ggtitle("") +
  theme(plot.title = element_text(size = 20)) +
  theme(axis.text=element_text(size=10)) + theme_classic() + theme(legend.text = element_text(size = 18))

# Define number of predictions to make
# Read the new data with deseasonalized and true volatility estimates 
data_full<-fread("eur_gbp_data_with_volatility.csv")

nobs <- dim(data_full)[1]

S_t_r <- as.numeric(model_return_robust$coefficients[1]) * data_full$D1 +
  as.numeric(model_return_robust$coefficients[2]) * data_full$D2 +
  as.numeric(model_return_robust$coefficients[3]) * data_full$D3 +
  as.numeric(model_return_robust$coefficients[4]) * data_full$D4 +
  as.numeric(model_return_robust$coefficients[5]) * data_full$D5 +
  as.numeric(model_return_robust$coefficients[6]) * data_full$D6 +
  as.numeric(model_return_robust$coefficients[7]) * data_full$D7 +
  as.numeric(model_return_robust$coefficients[8]) * data_full$D8 +
  as.numeric(model_return_robust$coefficients[9]) * data_full$D9 +
  as.numeric(model_return_robust$coefficients[10]) * data_full$D10 +
  as.numeric(model_return_robust$coefficients[11]) * data_full$D11 +
  as.numeric(model_return_robust$coefficients[12]) * data_full$D12 +
  as.numeric(model_return_robust$coefficients[13]) * data_full$D13 +
  as.numeric(model_return_robust$coefficients[14]) * data_full$D14 +
  as.numeric(model_return_robust$coefficients[15]) * data_full$D15 +
  as.numeric(model_return_robust$coefficients[16]) * data_full$D16

ReturnDeseasonalized <- data_full$ReturnStandardized - S_t_r

true_epsilon <- c(arima12$residuals, rep(NA, nobs - n_train))

for (i in (n_train+1):nobs) {
  true_epsilon[i] <- ReturnDeseasonalized[i] - arima12$coef["ar1"] * ReturnDeseasonalized[i-1] - arima12$coef["ar2"] * ReturnDeseasonalized[i-2] - arima12$coef["ma1"] * true_epsilon[i-1]
}

log_S_epsilon_squared <- as.numeric(model_volatility_robust$coefficients[1]) * data_full$D1 +
  as.numeric(model_volatility_robust$coefficients[2]) * data_full$D2 +
  as.numeric(model_volatility_robust$coefficients[3]) * data_full$D3 +
  as.numeric(model_volatility_robust$coefficients[4]) * data_full$D4 +
  as.numeric(model_volatility_robust$coefficients[5]) * data_full$D5 +
  as.numeric(model_volatility_robust$coefficients[6]) * data_full$D6 +
  as.numeric(model_volatility_robust$coefficients[7]) * data_full$D7 +
  as.numeric(model_volatility_robust$coefficients[8]) * data_full$D8 +
  as.numeric(model_volatility_robust$coefficients[9]) * data_full$D9 +
  as.numeric(model_volatility_robust$coefficients[10]) * data_full$D10 +
  as.numeric(model_volatility_robust$coefficients[11]) * data_full$D11 +
  as.numeric(model_volatility_robust$coefficients[12]) * data_full$D12

u_hat <- log(true_epsilon^2) - log_S_epsilon_squared

S_epsilon_squared <- rep(NA, nobs)

for (t in 1:nobs) {
  S_epsilon_squared[t] <- exp(log_S_epsilon_squared[t]) * sum(exp(u_hat[1:t])) / t
}

epsilon_star1 <- c(epsilon_star, rep(NA, nobs-n_train))

for (t in (n_train+1):nobs) {
  epsilon_star1[t] <- true_epsilon[t]/sqrt(S_epsilon_squared[t]) - coef(m2)["mu"]
}

sigma_star_squared <- rep(NA, nobs)
sigma_star_squared[n_train]<-var(data_full$ReturnStandardized[1:n_train])

for (i in (n_train+1):nobs) {
  sigma_star_squared[i] <- coef(m2)["omega"] + coef(m2)["alpha1"] * (epsilon_star1[i-1])^2 + coef(m2)["beta1"] * sigma_star_squared[i-1]
}

deseasonalized_volatility <- c(volatility(m2), sqrt(sigma_star_squared[(n_train+1):nobs]))

volatility_estimated <- c(volatility(m2), rep(NA, length(sqrt(sigma_star_squared[(n_train+1):nobs]))))
volatility_forecasted <- c(rep(NA, length(volatility(m2))), sqrt(sigma_star_squared[(n_train+1):nobs]))

# Compute unconditional volatility
volatility_unconditional <- sqrt( coef(m2)["omega"] / (1 - coef(m2)["alpha1"] - coef(m2)["beta1"]) )
volatility_unconditional <- rep(volatility_unconditional, nobs)
index <- seq(1, nobs, 1)

# Create a dataframe with estimated, forecasted and unconditional deseasonalized volatility
df_volatility_to_plot <- as.data.frame(cbind(index, volatility_estimated, volatility_forecasted, volatility_unconditional))

# Plot Deseasonalized Volatility of Standardized EUR/GBP 5-minute Log Returns
ggplot(df_volatility_to_plot, aes(x=index)) +
  labs(y = "", x = "") +
  geom_line(aes(y=volatility_estimated, col = "Estimated deseasonalized volatility", linetype = "Estimated deseasonalized volatility")) +
  geom_line(aes(y=volatility_forecasted, col = "Forecasted deseasonalized volatility", linetype = "Forecasted deseasonalized volatility")) +
  geom_line(aes(y=volatility_unconditional, col = "Unconditional deseasonalized volatility", linetype = "Unconditional deseasonalized volatility")) +
  scale_color_manual(name = "", breaks = c("Estimated deseasonalized volatility", "Forecasted deseasonalized volatility", "Unconditional deseasonalized volatility"), values = c("steelblue", "darkgreen", "darkred")) +
  scale_linetype_manual(name = "", breaks = c("Estimated deseasonalized volatility", "Forecasted deseasonalized volatility", "Unconditional deseasonalized volatility"), values = c("solid", "solid", "solid")) +
  ggtitle("") +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(axis.text=element_text(size=10)) + theme_classic() + theme(legend.text = element_text(size = 18)) + ylim(0,5)

true_volatility_estimated <- volatility(m2)*sqrt(S_epsilon_squared[1:n_train])
true_volatility_forecasted <- sqrt(sigma_star_squared[(n_train+1):nobs]) * sqrt(S_epsilon_squared[(n_train+1):nobs])

true_volatility <- c(true_volatility_estimated, true_volatility_forecasted)

true_volatility_estimated1 <- c(true_volatility_estimated, rep(NA, length(true_volatility_forecasted)))
true_volatility_forecasted1 <- c(rep(NA, length(true_volatility_estimated)), true_volatility_forecasted)

# Create a dataframe with estimated, forecasted and unconditional deseasonalized volatility
df_true_volatility_to_plot <- as.data.frame(cbind(index, true_volatility_estimated1, true_volatility_forecasted1))

# Plot True Volatility of Standardized EUR/GBP 5-minute Log Returns
ggplot(df_true_volatility_to_plot, aes(x=index)) +
  labs(y = "", x = "") +
  geom_line(aes(y=true_volatility_estimated1, col = "Estimated true volatility", linetype = "Estimated true volatility")) +
  geom_line(aes(y=true_volatility_forecasted1, col = "Forecasted true volatility", linetype = "Forecasted true volatility")) +
  scale_color_manual(name = "", breaks = c("Estimated true volatility", "Forecasted true volatility"), values = c("steelblue", "darkgreen")) +
  scale_linetype_manual(name = "", breaks = c("Estimated true volatility", "Forecasted true volatility"), values = c("solid", "solid")) +
  ggtitle("") +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(axis.text=element_text(size=10)) + theme_classic() + theme(legend.text = element_text(size = 18)) + ylim(0,10)

# Save volatility estimates and forecasts
volatility_final <- as.data.frame(cbind(deseasonalized_volatility, true_volatility, rep(volatility_unconditional[1], length(true_volatility)), seq(1, length(true_volatility), 1)))
