exprs2fingerprint <-
function(exprs, platform, species, progressBar = TRUE)
# Produces a normalized fingerprint directly from an exprs object
# GEOnormalize is used to normalize the score against the GEO background
# The pathways in the fingerprint are a combination of 
# * KEGG
# * Wikipathways
# * KEGG shell
# * Netpath
# * Static networks
#
# Platforms should be of the type listed in GEO (e.g. "GPL570").
# Species can be "mouse" or "human"
# Updated for v0.3
  {

#######
# Load data and check parameters
#######

# Load annotation and background data associated with each platform   
#    if (!(exists("chipframe"))){
#      print("Loading background file")
#      data(chipframe)
#    }
# check platform is compatible
	if (!(platform %in% names(chipframe))) stop("Platform name is invalid or not currently supported")

# Load pathway gene sets and check species

#	 if (!(exists("genesets")))
#	 { 	
#      data(genesets)
#    }
	if (!(species %in% names(genesets))) stop("Species name invalid or not supported")
#	data(list = genesets[species])
	gsdb <- get(genesets[species]) 
	
	
#########							
# Run pathway enrichment	
#########

    print("Running fingerprint")

# Annotate array according to the annotations
    geo.ann <- customCDFAnn(exprs, chipframe[[platform]]$ann)
    
# Use SCE to calculate a score for each pathway based on the mean or median
		
	SCE <- single.chip.enrichment( 
   			exprs = geo.ann,
  			geneset = gsdb,
			transformation = "squared.rank",
			statistic = "mean",
			normalizedScore = FALSE,
			progressBar = progressBar
   			)

#########							
# Threshold according to GEO corpus background 	
#########   			

	SCE.threshold<-thresholdFingerprint(SCE = SCE, platform = platform)
	
    return(SCE.threshold)

# End exprs2fingerprint
  }

