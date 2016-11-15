gsource("c:/CGP2009/Oliver/OPAM.library.R")     # load LIBRARY

# Project medulloblastoma outcome dataset from Pomeroy et al 2002 (dataset "C") in the space of 
# a few selected outcome-related gene sets from MSigDB release 2.5

outcome.gene.sets <- c("YU_CMYC_UP", "MYC_ONCOGENIC_SIGNATURE", "MENSSEN_MYC_UP", "METASTASIS_ADENOCARC_UP", 
                       "POMEROY_MD_TREATMENT_GOOD_VS_POOR_DN", "CHANG_SERUM_RESPONSE_UP")
OPAM.project.dataset(
   input.ds =                "c:/CGP2009/Oliver/Dataset_C_MD.gct",
   output.ds =               "c:/CGP2009/Oliver/Dataset_C_MD.pathways.gct",
   gene.set.databases =      c("c:/CGP2009/Oliver/MSIGDB_c2.all.v2.5.symbols.gmt"),
   gene.set.selection =      outcome.gene.sets,
   sample.norm.type =        "rank",  # "rank", "log" or "log.rank"
   output.score.type =       "ES")  # "ES" or "NES"

# Now sort projected dataset according to ROC with Relapse vs No-relapse phenotype 
# (ROC and p-value on right hand side of heatmap)

OPAM.match.projection.to.phenotypes(
    input.ds =               "c:/CGP2009/Oliver/Dataset_C_MD.pathways.gct",
    input.cls =              "c:/CGP2009/Oliver/Dataset_C_MD.cls",
    results.dir =            "c:/CGP2009/Oliver/",
    normalize.score =        T,
    normalization.type =     "zero.one",
    markers.num =            20,
    user.colors =            NA,
    markers.metric =         "ROC",   # "ROC" or "T.TEST"
    markers.file =           "c:/CGP2009/Oliver/Dataset_C_MD.markers.gct",
    sort.phenotypes =        T,
    sort.decreasing =        F,    # T = decreasing, F = increasing
    legend =                 T,
    char.res =               0.9,
    only.up =                T)

# Project dataset in the space of the entire database

OPAM.project.dataset(
   input.ds =                "c:/CGP2009/Oliver/Dataset_C_MD.gct",
   output.ds =               "c:/CGP2009/Oliver/Dataset_C_MD.pathways.ALL.gct",
   gene.set.databases =      c("c:/CGP2009/Oliver/MSIGDB_c2.all.v2.5.symbols.gmt"),
   gene.set.selection =      "ALL",
   sample.norm.type =        "rank",  # "rank", "log" or "log.rank"
   output.score.type =       "ES")  # "ES" or "NES"

# Sort projected dataset according to ROC with Relapse vs No-relapse phenotype 
# (ROC and p-value on right hand side of heatmap)

OPAM.match.projection.to.phenotypes(
    input.ds =               "c:/CGP2009/Oliver/Dataset_C_MD.pathways.ALL.gct",
    input.cls =              "c:/CGP2009/Oliver/Dataset_C_MD.cls",
    results.dir =            "c:/CGP2009/Oliver/",
    normalize.score =        T,
    normalization.type =     "zero.one",
    markers.num =            10,    # select the top 10/bottom 10 pathways
    user.colors =            NA,
    markers.metric =         "ROC",   # "ROC" or "T.TEST"
    markers.file =           "c:/CGP2009/Oliver/Dataset_C_MD.markers.gct",
    sort.phenotypes =        T,
    sort.decreasing =        F,    # T = decreasing, F = increasing
    legend =                 T,
    char.res =               1.25,
    only.up =                F)


