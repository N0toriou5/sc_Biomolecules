

### Area: Analytical Rank EnrichmeNt Analysis
area<-function(
    signatures,
    groups,
    sweights=NULL,
    gweights=NULL,
    minsize=1
){
    ### Convert to list if groups are not a list
    if(!is.list(groups)){
        groups<-apply(groups,1,function(x){
            hits<-colnames(groups)[x==1]
            return(hits)
        })
    }
    
    ### Remove from groups what is not in signature matrix
    features<-rownames(signatures)
    groups<-lapply(groups,function(x){
        y<-intersect(x,features)
        return(y)
    })
    
    ### Remove small groups
    groups<-groups[sapply(groups,length)>=minsize]
    
    ### Treat single "signature"
    if (is.null(nrow(signatures))){
        signatures <- matrix(signatures, length(signatures), 1, dimnames=list(names(signatures), "sample1"))
    } 
    
    
    ### Generate dummy signature weights
    if(is.null(sweights)){
        sweights<-matrix(1,nrow=nrow(signatures),ncol=ncol(signatures))
        dimnames(sweights)<-dimnames(signatures)
    }
    if(!identical(dim(signatures), dim(sweights))){
        stop("Signatures and Signature weights must be matrices of identical size")
    }
    
    ### Generate dummy group weights
    if(is.null(gweights)){
        gweights<-relist(rep(1,sum(sapply(groups,length))),skeleton=groups)
    }
    
    ### Apply weights to group belonging
    wgroups<-gweights
    for(i in 1:length(wgroups)){
        names(wgroups[[i]])<-groups[[i]]
    }
    rm(gweights)
    
    ### Rank-transform columns
    ranks<-apply(signatures,2,rank,na.last="keep")
    ### Assign a 0 to signature weights where the signature was NA
    sweights[is.na(ranks)]<-0
    ### 0-1 bound ranks
    boundranks<-t(t(ranks)/(colSums(!is.na(signatures))+1))
    ### Treat bound ranks as quantiles in a gaussian distribution (0=-Inf, 1=+Inf)
    gaussian <- qnorm(boundranks)
    ### Deal with NAs
    gaussian[is.na(gaussian)]<-0
    
    ### Apply signature weights to the normalized distribution
    gaussian<-gaussian*sweights
    
    
    ### Next, we see how each of the groups are behaving in these normalized signatures
    ### Create a boolean matrix with ngroup columns and signaturelength rows, indicating the matches 
    matches <- sapply(wgroups, function(group, allElements) {
        hereMatches<-as.integer(allElements%in%names(group))
        names(hereMatches)<-allElements
        # Weigth by group belonging
        weightedMatches<-hereMatches
        weightedMatches[names(group)]<-weightedMatches[names(group)]*group
        return(weightedMatches)
    }, allElements=rownames(gaussian))
    # And then transpose it
    matches<-t(matches)
    colnames(matches)<-rownames(signatures)
    
    # Number of matches per group
    groupmatches <- rowSums(matches)
    
    # Relative part of the signature that matches
    relativematches<-matches/groupmatches
    
    # This trick will overweight massively small groups with all their components highly-ranked.
    # Extreme case is with a group with one gene at the top
    
    # The core linear algebra operation. The true magic of rea
    enrichmentScore <- relativematches %*% gaussian
    
    # Finally, every enrichment is square-rooted to respect the criterion of normality
    normalizedEnrichmentScore<-enrichmentScore*sqrt(groupmatches)
    
    # Return output
    return(normalizedEnrichmentScore)
}



### Test it
if(FALSE){
    ### Signature matrix and Groups to match to the matrix
    set.seed(1)
    signatures<-matrix(runif(1000),ncol=10)
    signatures[,1]<-100:1/100
    colnames(signatures)<-LETTERS[1:10]
    rownames(signatures)<-paste0("row",1:nrow(signatures))
    groups<-list(GROUP1=c("row1","row2","row3"),GROUP2=c("row10","row2","row34","row4"),GROUP3="row1")
    
    # With some extra NAs
    signatures[sample(length(signatures),100)]<-NA
    
    # With group belonging weight
    gweights<-relist(runif(sum(sapply(groups,length))),skeleton=groups)
    
    # With signature weight
    sweights<-matrix(runif(length(signatures)),nrow=nrow(signatures),ncol=ncol(signatures))
    
    
    
    ### Run arena
    nes<-area(signatures=signatures,groups=groups,sweights=sweights,gweights=gweights,minsize=1)
    
    
    
    
}



