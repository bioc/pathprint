# Example script for hyperPathway

library(GMAfunctions)

# take gene sets from pathprint package
library(pathprint)
data(pathprint.Hs.gs)


# test with a random genelist of Entrez Gene IDs
sampleGenes<-c("5692", "5682", "29127", "5702", "84950", "1659", "5684", "5431", 
"115106", "11157", "24148", "23165", "11269", "22938", "5708", 
"10213", "10403", "79680", "5436", "57819", "9861", "55142", 
"79441", "91746", "58509", "988", "7251", "5422", "59286", "890", 
"83540", "11316", "5683", "5713", "5901", "57461", "5515", "5700", 
"55696", "85002", "51639", "151903", "26036", "23423", "25799", 
"6319", "84439", "1213", "2597", "57620", "55244", "699", "199704", 
"55735", "64965", "57835", "64852", "57465", "127703", "10713", 
"51691", "388531", "154007", "1315", "51477", "84922", "64236", 
"29882", "4247", "51264", "3805", "23474", "60314", "26994", 
"150946", "4738", "57474", "147015", "3190", "9055", "57642", 
"5435", "11325", "64763", "55258", "90407", "124446", "27250", 
"149175", "9821", "91869", "253635", "65243", "286826", "55083", 
"84975", "79622", "6457", "9219", "5128", "29945", "27316", "56931", 
"253827", "51473", "55957", "1733", "821", "79818", "115677", 
"5252", "115361", "8291", "23064", "5146", "837", "2779", "3842", 
"57456", "81853", "57419", "255104", "57117", "9787", "286527", 
"3763", "8621", "27121", "1009", "4832", "89777", "4504", "5961", 
"4170", "2155", "91978", "57082", "192134", "154075", "3274", 
"4706", "83938", "147945", "2838", "26115", "25900", "79443", 
"10087", "121441", "54039", "10594", "10114", "5211", "2785")

########
# if have gene symbols rather than entrez IDs
# can use bioconductor annotation package to convert to Entrez Gene IDs
symbols<-c("PSMB4", "PSMA1", "RACGAP1", "PSMC3", "PRPF38A", 
"DHX8", "PSMA3", "POLR2B", "HAUS1", "LSM6", "PRPF6", "NUP205", 
"DDX19B", "SNW1", "PSMD2", "PSMD14", "NDC80", "C22orf29", "POLR2G", 
"LSM2", "PSMD6", "HAUS2", "HAUS3", "YTHDC1", "C19orf29", "CDC5L", 
"TSG101", "POLA1", "UBL5", "CCNA2", "NUF2", "COPE", "PSMA2", 
"PSMD7", "RAN", "ISY1", "PPP2CA", "PSMC1", "RBM22", "FAM86B1", 
"SF3B14", "CCDC12", "ZNF451", "TMED3", "ZNF324", "SCD", "HHIPL1", 
"CLTC", "GAPDH", "STIM2", "SLC47A1", "BUB1", "ZNF585A", "DNAJC11", 
"MRPS9", "SLC4A5", "TUT1", "TBC1D24", "C1orf216", "USP39", "NAA38", 
"RGS9BP", "SNRNP48", "COPB1", "ISYNA1", "FIZ1", "PDLIM2", "ANAPC2", 
"MGAT2", "MRPL27", "KIR2DL4", "ETHE1", "C12orf10", "RNF11", "FAM59B", 
"NEDD8", "ZNF490", "DHRS13", "HNRNPK", "PRC1", "COL20A1", "POLR2F", 
"DDX42", "ZNF574", "THNSL2", "TMEM41A", "TMEM219", "PDCD4", "MANEAL", 
"RB1CC1", "RFT1", "CCDC75", "ZNF643", "LIN9", "KIF26B", "MFSD5", 
"SNRNP25", "SH3GL3", "MTA2", "CDK17", "ANAPC4", "RBMX", "DUS3L", 
"MSRB3", "DCDC2", "LIN37", "DIO1", "CANX", "ZNF552", "NOSTRIN", 
"PHF1", "GBP4", "DYSF", "SETX", "PDE6C", "CASP4", "GNAT1", "TNPO1", 
"KIAA1143", "TMEM14B", "SLC24A3", "TMCO4", "INTS12", "DLGAP5", 
"TMSB15B", "KCNJ6", "CDK13", "DKK4", "CDH11", "NME3", "SERPINB12", 
"MT3", "PRPH2", "MCL1", "F7", "C19orf20", "CASC5", "B3GNT6", 
"SAMD3", "HRH2", "NDUFAB1", "C10orf11", "NLRP4", "GPR15", "TANC2", 
"IFFO1", "FYCO1", "COL4A3BP", "NEDD1", "PCBP3", "PRPF8", "HIPK3", 
"PFKL", "GNG3")

library(org.Hs.eg.db)
entrez.gene.ids<-unlist(mget(symbols, org.Hs.egSYMBOL2EG, ifnotfound = NA))
##########

# Run pathway enrichment
# estimate approximate genenome-wide background, 20000 genes
pathwayEnrichment<-hyperPathway(genelist = sampleGenes,
								geneset = pathprint.Hs.gs,
								Nchip = 20000
								)

# examine top 10 by p-value
# (don't view pathwayname, column 5, as same as rownames)

head(pathwayEnrichment[order(pathwayEnrichment$"P-value", decreasing = FALSE),-6],10)

################
# pathwayHeatmap example, either using a p-value cutoff
# function to create a heatmap of the enriched pathways, either cutoff by p-value or numnber of genes
####
# List of arguments
# genelist - vector of genes
# geneset - list of pathways or genesets
# cutoff - "P-value" or "nGenes", method by which to determine which pathways to show
# pthreshold - P-value cutoff
# gthreshold - nGenes cutoff
# convertEntrez2Symbol - whether or not to convert Entrez gene ids to symbols for annotation
# species - species for conversion
#####
# Example using a pathway cutoff
pathway.matrix<-pathwayHeatmap(genelist = sampleGenes,
							   geneset = pathprint.Hs.gs,
							   cutoff = "P-value",
							   pthreshold = 0.001,
							   Nchip = 20000,
							   convertEntrez2Symbol = TRUE,
                         	   species = "human",
                         	   mar = c(10,15), # further arguments are passed onto heatmap
                         	   cexRow = 0.5, cexCol = 0.5
                         	   )
# or a number of genes cutoff
pathway.matrix<-pathwayHeatmap(genelist = sampleGenes,
							   geneset = pathprint.v0.3.Hs.gs,
							   cutoff = "nGenes",
							   gthreshold = 10,
							   Nchip = 20000,
							   convertEntrez2Symbol = TRUE,
                         	   species = "human",
                         	   mar = c(10,15), # further arguments are passed onto heatmap
                         	   cexRow = 0.5, cexCol = 0.5
                         	   )
