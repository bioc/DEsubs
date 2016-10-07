.constructNetwork <- function( org, mRNAexpr, mRNAnomenclature='entrezgene', 
                                pathways='all')
{
    # If a filename is given, read a text file from the 'User' directory.
    if ( class(mRNAexpr) == 'character' )
    {
        dir <- paste0(cache[['baseDir']], '//User') 
        mRNAexpr <- readLines(paste(dir, mRNAexpr, sep='//'))
        mRNAexpr <- strsplit(mRNAexpr, '\t')
        mRNAexpr <- do.call('rbind', mRNAexpr)
        rownames(mRNAexpr) <- mRNAexpr[, 1]
        mRNAexpr <- mRNAexpr[, -1]
        class(mRNAexpr) <- 'numeric' 
    }

    # Change nomenclature if necessary

    if ( mRNAnomenclature != 'entrezgene' ) 
    {
        file <- 'extdata//Data//libraryEntrezToExternalNomenclature.RData'
        file  <- system.file(file, package='DEsubs')
        load(file, e <- new.env())
        orgLib <- e[['libraryEntrezToExternalNomenclature']][[org]]
        lib <- orgLib[[mRNAnomenclature]]

        # Keep rows with entrez ids
        idx <- which( as.character(rownames(mRNAexpr)) %in% lib[, 2] )
        mRNAexpr <- mRNAexpr[idx, ]
        
        # Convert to entrez ids
        xlib <- lib[, 1]
        names(xlib) <- lib[, 2]
        rownames(mRNAexpr) <- xlib[rownames(mRNAexpr)]        
    }

    # Create pathway graph
    file  <- system.file('extdata//Data//edgeLists.RData', package='DEsubs')
    load(file, e <- new.env())
    edgeList <- e[['edgeLists']][[org]]

    # Filter edgelist according to pathway type
    idx <- sapply(edgeList, is.factor)
    edgeList[idx] <- lapply(edgeList[idx], as.character)

    pathIds <- as.numeric(edgeList[, 4])
    if ( pathways == 'Metabolic')
        { edgeList <- edgeList[which(pathIds < 2000), ] }
    if ( pathways == 'Non-Metabolic')
        { edgeList <- edgeList[which(pathIds >= 2000), ] }

    # Filter edgelist with gene within the RNA-seq data
    idx1 <- which(edgeList[, 1] %in% rownames(mRNAexpr))
    idx2 <- which(edgeList[, 2] %in% rownames(mRNAexpr))
    idx <- intersect(idx1, idx2)
    edgeList <- edgeList[idx, ]


    return( list(mRNAexpr=mRNAexpr, edgeList=edgeList) )
}


.pruneNetwork <- function(   edgeList, mRNAexpr, DEGchoice, classes, 
                            DEGthresh, corr_threshold , org, verbose=TRUE, 
                            CORtool, rankedList=NULL)
{

    CORtool.options <- c('pearson', 'kendall', 'spearman')
    if ( !CORtool %in% CORtool.options)
    {
        message(CORtool, ' option not availiable.')
        message('Availiable options:', CORtool.options)
    }

    #
    # Phase 1 - Prune nodes
    #
    if ( verbose ) { message('Pruning nodes...', appendLF = FALSE) }

    # Apply a built-in differential expression analysis method to extract 
    # a list of differentially expressed genes (Q-values).
    if ( is.null(rankedList) )
    {
        KEGGgenes <- unique(as.vector(t(edgeList[, 1:2])))
        keggIdx <- which( rownames(mRNAexpr) %in% KEGGgenes )        
        mRNAexpr <- mRNAexpr[keggIdx, ]
        genes <- .DEanalysis( count.matrix=mRNAexpr, DEGchoice=DEGchoice, 
                            classes=classes )
    }
    # The Q-values for the differentially expressed genes have been supplied
    # as input and no buit-in method is used, even if it selected
    if ( !is.null(rankedList) )
    {
        # A filepath is given as input
        if ( class(rankedList) != 'numeric')
        {
            file <- rankedList
            rankedList.data <- read.table(file, dec = '.', sep=' ')
            genes <- as.numeric(rankedList.data[, 2]) # Q-values
            names(genes) <- rankedList.data[, 1] # Gene names
        }
        # A ranked list of differentially expressed genes is given as input,
        # in the form of a named vector, storing the Q-values and the gene
        # names (the vector's names).
        if ( class(rankedList) == 'numeric')
        { 
            genes <- rankedList 
        }
    }

    # if (verbose)
    #     { message('\n\tGenes before NodeScore: ', length(genes)) }
    
    # Keep statistically significant degs
    degs <- genes[genes < DEGthresh]
    
    # if (verbose)
    #     { message('\tGenes after NodeScore: ', length(degs)) }

    # if (verbose)
    #     { message('\tEdges before NodeScore: ', nrow(edgeList)) }


    # Filter edgelist with degs by discarding all non significant degs.
    idx1 <- which(edgeList[, 1] %in% names(degs)) # Both nodes of an edge
    idx2 <- which(edgeList[, 2] %in% names(degs)) # have to be degs
    edgeList <- edgeList[intersect(idx1, idx2), ] 

    # if (verbose)
    #     { message('\tEdges after NodeScore: ', nrow(edgeList)) }

    
    if ( verbose ) { message('done.', appendLF = TRUE) }

    #
    # Phase 2 - Prune edges
    #

    if ( verbose ) 
        { message('Pruning edges...', appendLF = FALSE) }
    # if (verbose)
    #     { message('\n\tEdges before EdgeScore: ', nrow(edgeList)) }
    
    # Prune edgelist
    edgeList <- .pruneEdges(edgeList, mRNAexpr, corr_threshold, CORtool)
    
    # Remove degs whose edges have been competely removed
    uGenes <- unique(as.vector(as.matrix(edgeList[, -3])))
    degs <- degs[ names(degs) %in% uGenes ]

    
    # if ( verbose )
    #     { message('\tEdges after EdgeScore: ', nrow(edgeList)) }
    if ( verbose ) 
        { message('done.', appendLF = TRUE) }

    return( list( 'edgeList'=edgeList, 'DEgenes'=degs) )
}


.pruneEdges   <- function( edgeList, mRNAexpr, thr, CORtool )
{
    # Find expression
    mapper <- seq_len(nrow(mRNAexpr))
    names(mapper) <- rownames(mRNAexpr)

    exprEdgelistSource <- mRNAexpr[mapper[edgeList[, 1]],]
    exprEdgelistDestin <- mRNAexpr[mapper[edgeList[, 2]],]

    # No edges
    if (is.null(edgeList) || nrow(edgeList) == 0)
    {
        return(matrix(, nrow=0, ncol=3))
    }

    CR <- vector(mode='numeric', length=nrow(edgeList))
    for ( i in seq_len(nrow(edgeList)) )
    {
        CR[i] <- cor(exprEdgelistSource[i,], exprEdgelistDestin[i,], 
                    method = CORtool )
    }

    # Add one extra column to edgelist for correlation score
    edgeList <- cbind(edgeList, 'corr'=CR)

    # Score edges
    hitIdx <- vector(mode='numeric', length=nrow(edgeList))

    for ( i in seq_len(nrow(edgeList)) )
    {
        type <- edgeList[i, 3]
        if ( type == 1 ) { reg <- 1 }
        if ( type == 2 ) { reg <- -1 }          
        if ( type != 3 && CR[i] * reg > thr ) { hitIdx[i] <- 1 }
    }

    # Keep edges passing the previous criteria 
    edgeListPruned <- edgeList[which(hitIdx == 1), , drop=FALSE]
    
    # Remove same edges coming from different pathways
    idx <- which(!duplicated(edgeListPruned[, c(1,2)]))
    edgeListPruned <- edgeListPruned[idx, ]

    if ( nrow(edgeListPruned) == 0 )
        { return( edgeListPruned ) }
    rownames(edgeListPruned) <- seq_len(nrow(edgeListPruned))
    
    
    return(edgeListPruned)
}

