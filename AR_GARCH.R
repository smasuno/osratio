wd<-"~/Desktop/FIN580 Project/"
setwd(wd)
require(forecast)
require(dynlm)
#require(quantmod)
require(vars)
require(coefplot)
require(data.table)
require(rugarch)

filenames <- list.files("~/Desktop/FIN580 Project/trading_strategy_data/", pattern="*.csv", full.names=TRUE)
coef_mat=data.frame(OBC_S.l1=numeric(),OSC_S.l1=numeric(),CBC_S.l1=numeric(),CSC_S.l1=numeric(),OBP_S.l1=numeric(),
                    OSP_S.l1=numeric(),CBP_S.ls=numeric(),CSP_S.l1=numeric(),VIX_DIFF.l1=numeric(),RET.l1=numeric(),
                    const=numeric())
ticker_vec<-vector()
counter=0;

#Getting Weekly Volumes
weekly_vol<-read.csv('Part4_weekly.csv',header=T)
weekly.table<-data.table(weekly_vol)
average_weekly_volume<-weekly.table[!is.na(NUMTRD),]
average_weekly_volume<-average_weekly_volume[,mean(NUMTRD),by='TICKER']

for(i in 1:length(filenames)){
temp<-read.csv(filenames[i])
vix<-read.csv('~/Desktop/FIN580 Project/VIX.csv',header=T)
vix<-vix[match(as.Date(as.character(temp$DATE),'%Y%m%d'),as.Date(vix$Date)),]
vix_ret<-diff(log(vix$Adj.Close))
vix_ret[is.na(vix_ret)]<-0
#Squared VIX
vix_ret<-vix_ret*vix_ret
ticker<-temp$TICKER[1]

#Take out NA in NUMTRD
temp<-temp[!is.na(temp$NUMTRD),]
temp<-temp[order(-temp$DATE),]

if(dim(temp)[1] > 11){
  counter=counter+1;
#OS Ratios
OBC_S<-temp$OBC/temp$NUMTRD
OSC_S<-temp$OSC/temp$NUMTRD
CBC_S<-temp$CBC/temp$NUMTRD
CSC_S<-temp$CSC/temp$NUMTRD
OBP_S<-temp$OBP/temp$NUMTRD
OSP_S<-temp$OSP/temp$NUMTRD
CBP_S<-temp$CBP/temp$NUMTRD
CSP_S<-temp$CSP/temp$NUMTRD

data<-cbind(temp$DATE,OBC_S,OSC_S,CBC_S,CSC_S,OBP_S,OSP_S,CBP_S,CSP_S,temp$RET)
data<-data[-1,]
data<-cbind(data,vix_ret)
colnames(data)<-c('DATE','OBC_S','OSC_S','CBC_S','CSC_S','OBP_S','OSP_S','CBP_S','CSP_S','RET','VIX_DIFF')
data<-data.frame(data)
data<-data[!is.nan(data$OBC_S),]
data<-data[is.finite(data$OBC_S),]
n_obs<-dim(data)[1]

#Multivariate TS Modelling
#model_ts<-ts(data[c('OBC_S','OSC_S','CBC_S','CSC_S','OBP_S','OSP_S','CBP_S','CSP_S','VIX_DIFF','RET')],start=min(as.Date(temp$DATE)),end=max(as.Date(temp$DATE)))

OBC.l1<-ts(data$OBC_S)[2:n_obs]
OSC.l1<-ts(data$OSC_S)[2:n_obs]
CBC.l1<-ts(data$CBC_S)[2:n_obs]
CSC.l1<-ts(data$CSC_S)[2:n_obs]
OBP.l1<-ts(data$OBP_S)[2:n_obs]
OSP.l1<-ts(data$OSP_S)[2:n_obs]
CBP.l1<-ts(data$CBP_S)[2:n_obs]
CSP.l1<-ts(data$CSP_S)[2:n_obs]
VIX.l1<-ts(data$VIX_DIFF)[2:n_obs]
RET.l1<-ts(data$RET)[2:n_obs]
RET<-ts(data$RET)[1:(n_obs-1)]

#model<-VAR(model_ts)
#coef(model$varresult$RET)
#model<-dynlm(RET ~ L(OBC_S,1) + L(OSC_S,1) + L(CBC_S,1) + L(CSC_S,1) + L(OBP_S,1) + L(OSP_S,1) + L(CBP_S,1) + L(CSP_S,1) + L(VIX_DIFF,1) + L(RET,1),data=model_ts)
model<-lm(RET ~ OBC.l1 + OSC.l1 + CBC.l1 + CSC.l1 + OBP.l1 + OSP.l1 + CBP.l1 + CSP.l1 + VIX.l1 + RET.l1)
residuals<-model$residuals

do_garch=function(x){tryCatch({
#Model Residuals of AR(1) model as GARCH(1,1)
gfit<-garchFit(~garch(1,1),residuals,trace=FALSE)

#Simulate GARCH residuals
#garch_rand<-rnorm((n_obs-2),mean=0,sd=1)*sqrt(gfit@fit$series$h)

##Weighted Least Squares
#Long Run Variance Calculation from GARCH(1,1)
long_run_var<-gfit@fit$coef[2]/(1-gfit@fit$coef[3]-gfit@fit$coef[4])

if(is.finite(long_run_var)){
#Rescale Returns by Conditional Variance
wRET<-RET/sqrt(gfit@fit$series$h)
wOBC.l1<-OBC.l1/sqrt(gfit@fit$series$h)
wOSC.l1<-OSC.l1/sqrt(gfit@fit$series$h)
wCBC.l1<-CBC.l1/sqrt(gfit@fit$series$h)
wCSC.l1<-CSC.l1/sqrt(gfit@fit$series$h)
wOBP.l1<-OBP.l1/sqrt(gfit@fit$series$h)
wOSP.l1<-OSP.l1/sqrt(gfit@fit$series$h)
wCBP.l1<-CBP.l1/sqrt(gfit@fit$series$h)
wCSP.l1<-CSP.l1/sqrt(gfit@fit$series$h)
wVIX.l1<-VIX.l1/sqrt(gfit@fit$series$h)
wRET.l1<-RET.l1/sqrt(gfit@fit$series$h)
#Linear Regression on Rescale Returns
wlsmodel<-lm(wRET ~ wOBC.l1 + wOSC.l1 + wCBC.l1 + wCSC.l1 + wOBP.l1 + wOSP.l1 + wCBP.l1 + wCSP.l1 + wVIX.l1 + wRET.l1)
model<-wlsmodel

return(1)
}
},error=function(err){
  return(0)
})}

do_garch()
#Coefficient Plot for Return
path<-paste('/Users/Shintaro/Desktop/FIN580 Project/coeff_graphs/',ticker,'.png',sep='')
png(filename = path)
fig<-coefplot(model,title=ticker)
plot(fig)
dev.off()

#Storing in coefficient matrix
coefficients<-as.vector(coef(model))
ticker_vec[counter]<-as.character(ticker)
coef_mat[counter,]<-c(coefficients)

print(ticker)
}
else{
  
}
}



