RSC_thresholding = function(avgLoss,R){
  tgrid   <- seq(0.005, 0.995, by=0.01)  
  #avgLoss <- colMeans(FLOSSES)
  avgLoss <- avgLoss
  plot(tgrid, avgLoss, t='b')
  topt  <- tgrid[which.min(avgLoss)]
  RSC   <- R
  RSC[abs(RSC)<=topt] <- 0
  return(RSC)
}