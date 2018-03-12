customCDFAnn <-
function(data, ann)
# Annotates according to custom cdf, averaging to resolve redundancies
{
    # Assign probe ID column
    data<-as.data.frame(data)
    if ("ID_REF" %in% colnames(data)){
        print("Using ID_REF column instead of rownames as probe IDs")
        } else {
        data <- cbind(rownames(data), data)
        colnames(data)[1] <- "ID_REF"
        }
    data$ID_REF <- tolower(data$ID_REF)
        
    # Merge with annotation dataframe
    if (sum(c("ID", "AFFYID") %in% colnames(ann)) < 1) 
        stop ("annotation frame does not contain expected columns")
    
    if ("ID" %in% colnames(ann)) {
        ann$ID <- tolower(ann$ID)
        data.temp <- merge(ann, data, by.x="ID", by.y="ID_REF")
    }

    else if ("AFFYID" %in% colnames(ann)) {
        ann$AFFYID <- tolower(ann$AFFYID)
        data.temp <- merge(ann, data, by.x="AFFYID", by.y="ID_REF")
    }
    data.temp <- data.temp[,-1]

# convert to numeric as occasionally factors can complicate matters
    data.temp[,-1] <- apply(data.temp[,-1, drop = FALSE], 2, 
                            function(x){as.numeric(as.character(x))})

# aggregate Entrez IDs using mean
# edit Nov 27th 2011, changed to median as more appropriate for the 
# human ST chips.
# edit Nov 28th 2011, changed back to mean. Not necessary to make this 
# modification at this stage
    data.agg <- aggregate(data.temp[,-1], list(EntrezID=data.temp$EntrezID), 
                        mean)
    rownames(data.agg) <- sub("_at", "", data.agg$EntrezID)
    data.agg<-data.agg[,-1, drop = FALSE]
    return(data.agg)
# End of customCDFAnn
}

