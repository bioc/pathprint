consensusFingerprint <-
function(fingerprintframe, threshold)
  # Produce a consensus fingerprint. For each pathway, the consesnsus fingerprint is
  # 1 if either the mean of the fingerprints exceed a certain threshold,
  # -1 if they fall below -threshold,
  # 0 otherwise
  {
    if (threshold < 0) stop("threshold must be positive")

    mean.frame <- apply(fingerprintframe, 1, mean) # create mean fingerprint
    consensus <- rep(0, nrow(fingerprintframe)) # initiate consenssu vector

    if (!(length(consensus) == length(mean.frame))) stop("Error")

    names(consensus) <- names(mean.frame)
    consensus[mean.frame > threshold] <- 1 # fill with pathways passing threshold
    consensus[mean.frame < (-threshold)] <- (-1)

    return (consensus)
  }

