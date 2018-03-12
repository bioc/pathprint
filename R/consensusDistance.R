consensusDistance <-
function(consensus, fingerprintframe)
# function to calculate distance from a consensus fingerprint, accounting only 
# for significant pathawys in the consensus
# consensus - consensus fingerprint
# fingerprintframe - dataframe of fingerprints from which the distance will be 
# calculated
# returns vector with names = colnames(fingerprintFrame)
# scale by max distance (2 x length of the consensus) 
# returns a p-value, based on an assumption of a normally distributed 
# fingerprintframe
    {
        # subset so that only significant pathawys in the conensus are used
        consensus.sub <- consensus[abs(consensus) > 0]
        consensus.length <- length(consensus.sub)
        frame.sub <- fingerprintframe[abs(consensus) > 0,]
        # calculate the distance
        frame.distance <- frame.sub - consensus.sub
        # sum over pathways
        distance<-colSums(as.data.frame(abs(frame.distance)), na.rm = TRUE)
        # scale against max length
        print(paste("Scaling against max length,", (2 * consensus.length)))
        scaled.distance <- distance/(2 * consensus.length)
        # order and return
        sorted.scaled.distance<-sort(scaled.distance)
        # now calculate p-values, based on assumption that 
        # sorted.scaled.distance is normally distributed
        mean.dist<-mean(sorted.scaled.distance)
        sd.dist<-sd(sorted.scaled.distance)
        pvals<-pnorm(sorted.scaled.distance, mean = mean.dist, sd = sd.dist)
        return(data.frame(distance = sorted.scaled.distance, pvalue = pvals))
        }
