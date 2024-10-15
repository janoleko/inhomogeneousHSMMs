AIC_BIC = function(mod, n){
  if(!is.na(mod$minimum)){
    AIC = 2*(mod$minimum + length(mod$estimate))
    BIC = 2*mod$minimum + log(n)*length(mod$estimate)
  } else{
    AIC = 2*(mod$value + length(mod$par))
    BIC = 2*mod$value + log(n)*length(mod$par)
  }
  res = c(AIC, BIC)
  names(res) = c("AIC", "BIC")
  return(res)
}


get_dwell_times = function(states, k){
  ind = which(states == k)
  dwell = 0
  counter = 1
  for (i in 1:(length(ind)-1)){
    if(ind[i+1] == ind[i]+1){
      counter = counter+1
    } else{
      dwell = c(dwell, counter)
      counter = 1
    }
  }
  return(dwell[2:length(dwell)])
}