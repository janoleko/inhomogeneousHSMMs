
# function that can compute AIC and BIC based on outputs from nlm, optim and nlminb
AIC_BIC = function(mod, n){
  if(length(mod$minimum)>0){
    AIC = 2*(mod$minimum + length(mod$estimate))
    BIC = 2*mod$minimum + log(n)*length(mod$estimate)
  } else if(length(mod$value)>0){
    AIC = 2*(mod$value + length(mod$par))
    BIC = 2*mod$value + log(n)*length(mod$par)
  } else if(length(mod$objective)>0){
    AIC = 2*(mod$objective + length(mod$par))
    BIC = 2*mod$objective + log(n)*length(mod$par)
  }
  res = c(AIC, BIC)
  names(res) = c("AIC", "BIC")
  return(res)
}