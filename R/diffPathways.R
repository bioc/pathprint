diffPathways <-
function(fingerprints, fac, threshold)
# function to return pathawys consistently different between two groups 
# of fingerprints
# fingerprints is a matrix of fingerprints
# The number of columns in fingerprints should correspond to the length of fac
# fac is a vector of levels or characters corresponding to two groups over 
# which to conduct the comparison
    {
    # convert to factor in case not already in correct class
    fac<-as.factor(fac)
    
    # check 2 levels and fac length
    levels<-levels(fac)
    if (length(levels) != 2) stop("fac must contain 2 levels")
    if (length(fac) != ncol(fingerprints)) 
        stop ("fac length must equal nubmer of fingerprints")
    if ((threshold < 0 | threshold > 2)) 
        stop ("threshold value should be between 0 and 2")

    # group 1
    if (sum(fac == levels[1]) == 1){
        print("N.B. only 1 sample for group 1")
        mean.1<-fingerprints[,fac == levels[1]]
        }
    else {mean.1<-apply(fingerprints[,fac == levels[1]], 1, mean)}
    
    # group 2
    if (sum(fac == levels[2]) == 1){
        print("N.B. only 1 sample for group 2")
        mean.2<-fingerprints[,fac == levels[2]]
        }
    else {mean.2<-apply(fingerprints[,fac == levels[2]], 1, mean)}

    # calculate difference in means
    mean.diff<-mean.1-mean.2
    
    # return pathways with difference greater than threshold
    signifPathawys<-na.omit(names(mean.diff)[abs(mean.diff) > threshold])
    return(as.character(signifPathawys))
    }

