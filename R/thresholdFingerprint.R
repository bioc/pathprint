thresholdFingerprint <-
function(SCE, platform){
#  data(platform.thresholds)
  # apply threshold
  thresholdTable<-platform.thresholds[[platform]]
  SCE.threshold<-(SCE>thresholdTable[,"high"]) - (SCE<thresholdTable[,"low"])
  return(SCE.threshold)
  }

