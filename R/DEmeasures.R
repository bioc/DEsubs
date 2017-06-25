# Gene level measures
#

geneVisualization <- function ( DEsubs.out, 
                measures.topological='all', 
                measures.functional='all',
                measures.barplot=TRUE, 
                topGenes=10,
                colors.topological=colorRampPalette(c("white", "red"))(100),
                colors.functional=colorRampPalette(c("white", "red"))(100),
                colors.barplot='#C7EDFCFF',
                size.topological=c(5,4),
                size.functional=c(7,4),
                size.barplot=c(5,5),
                outfile.topological='',
                outfile.functional='',
                outfile.barplot='',
                export='pdf',
                verbose=TRUE)    
{
    if ( verbose )
        { message('Generating gene-level view...', appendLF = FALSE) }
    
    edgeList <- DEsubs.out[['edgeList']]
    DEgenes  <- DEsubs.out[['DEgenes']]
    org      <- DEsubs.out[['org']]
    nomen    <- DEsubs.out[['mRNAnomenclature']]

    if ( length(DEgenes) < topGenes ) { topGenes <- length(DEgenes) }

    # Keep top DEgenes
    topGenes <- sort(DEgenes)[1:topGenes]
    output   <- list()

    # Topological analysis
    supportedTopologicalMeasures <- c(  'degree', 
                                        'betweenness', 
                                        'closeness',
                                        'hub_score', 
                                        'eccentricity',
                                        'page_rank')
    # Ontological analysis
    supportedFunctionalMeasures <- .getExternalMeasures()

    if ( !is.null(measures.topological) )
    {
        if ( 'all' %in% measures.topological )
            { measures.topological <- supportedTopologicalMeasures }

        failTopo <- measures.topological[!measures.topological %in% 
                                        supportedTopologicalMeasures]

        if ( length(failTopo) > 0 )
            { message('Topological option(s) ', failTopo, ' not supported.') }
    }

    if ( !is.null(measures.functional) )
    {
        if ( 'all' %in% measures.functional )
            { measures.functional <- supportedFunctionalMeasures }

        failFunc <- measures.functional[!measures.functional %in% 
                                        supportedFunctionalMeasures]

        if ( length(failFunc) > 0 )
            { message('Ontological option(s) ', failFunc, ' not supported.') }
    }

    # Analysis #1
    topology <- .topologicalMeasures(measures=measures.topological,
                                    edgeList=edgeList, 
                                    topGenes=topGenes, 
                                    org=org,
                                    visualize=TRUE,
                                    colors=colors.topological,
                                    width=size.topological[1],
                                    height=size.topological[2],
                                    outfile=outfile.topological,
                                    export=export)

    # Analysis #2
    ontology <- .ontologicalMeasures(measures=measures.functional,
                                    topGenes=topGenes, 
                                    org=org,
                                    visualize=TRUE,
                                    colors=colors.functional,
                                    width=size.functional[1],
                                    height=size.functional[2],
                                    outfile=outfile.functional,
                                    export=export)
    # Analysis #3
    if (  measures.barplot )
    {
        if ( org == 'hsa' )
        {
            names(topGenes) <- .changeAnnotation(annData=names(topGenes), 
                                        org='hsa', choice='entrezToHGNC')
        }

        .matrixVisualization( -log10(topGenes), type='barplot', title='',
                            colors=colors.barplot, 
                            width=size.barplot[1],
                            height=size.barplot[2],
                            outfile=outfile.barplot,
                            export=export)
    }

    output <- c(output, list('measures.topological'=topology) )
    output <- c(output, list('measures.functional'=ontology) )

    if ( verbose )
        { message('done', appendLF = TRUE) }

    return( output )
}

.topologicalMeasures <- function( edgeList, measures, topGenes, org, 
                        visualize, colors, width, height, outfile, export)
{
    if ( is.null(measures) ) 
    { return(NULL) } 

    supportedMeasures <- c( 'degree', 
                            'betweenness', 
                            'closeness',
                            'hub_score', 
                            'eccentricity',
                            'page_rank')

    if ( missing (measures) ) { measures <- supportedMeasures }

    gi   <- graph_from_edgelist(as.matrix(edgeList[, 1:2]))
    topo <- list()
    topologicalHandlers <- .getTopologicalHandlers()
    for ( i in seq_len(length(measures)) )
    {
        res <- sort(topologicalHandlers[[measures[i]]](gi), decreasing=TRUE)
        topo <- c(topo, list(res))
    }
    names(topo) <- measures

    uGenes  <- unique(as.vector(as.matrix(edgeList[, 1:2])))
    ranking <- matrix(, nrow=length(uGenes), ncol=length(topo))
    for ( i in seq_len(length(topo)) )
    {
        ranking[, i] <- topo[[i]][as.character(uGenes)]
        ranking[, i] <- floor(ranking[, i]/max(ranking[, i])*100)/100
    }
    rownames(ranking) <- uGenes
    colnames(ranking) <- names(topo)


    # Return info about top genes
    if ( !missing(topGenes) && !is.null(topGenes) )
    {
        idx     <- which( rownames(ranking) %in% names(topGenes) )
        ranking <- ranking[idx, , drop=FALSE]
    }


    if ( visualize )
    {
        if (org == 'hsa')
        {
            rownames(ranking) <- .changeAnnotation( annData=rownames(ranking), 
                                                    org='hsa',
                                                    choice='entrezToHGNC')
        }

        .matrixVisualization( dataVis=as.matrix(ranking), 
                                type='heatmap', title='topological',
                                colors=colors,
                                width=width,
                                height=height,
                                outfile=outfile,
                                export=export )
    }

    return(ranking)
}

.ontologicalMeasures <- function( measures, topGenes, org, visualize, colors,
                        width, height, outfile, export )
{
    if ( org != 'hsa' ) 
        { return(NULL) } # Homo sapiens only
    if ( is.null(measures) ) 
        { return(NULL) } 

    supportedTerms <- .getExternalMeasures()

    if ( missing (measures) ) { measures <- supportedTerms }

    # Change gene annotation to HGNC for hsa 
    names(topGenes) <- .changeAnnotation(annData=names(topGenes), org='hsa',
                                        choice='entrezToHGNC')

    ranking <- matrix(, nrow=length(topGenes), ncol=length(measures))
    invalidIdx <- c()
    for ( j in seq_len(length(measures)) )
    {
        targets <- .loadTermData( type=measures[j] )
        if ( is.null(targets) )
        {
            invalidIdx <- c(invalidIdx, j)
            next()
        }
        for ( i in seq_len(length(topGenes)) )
        {
            idx <- targets %in% names(topGenes[i])
            idx <- matrix( idx, nrow=nrow(targets) )*1
            ranking[i, j] <- length(which(rowSums( idx ) > 0))
        }
        if ( max(ranking[, j]) > 0 )
        {
            ranking[, j] <- floor(ranking[, j]/max(ranking[, j])*100)/100
        }
    }
    rownames(ranking) <- names(topGenes)
    colnames(ranking) <- measures
    if ( length(invalidIdx) > 0 ) { ranking <- ranking[, -invalidIdx]}

    if ( visualize )
    {
        .matrixVisualization( dataVis=(as.matrix(ranking) ), 
                            type='heatmap', title='external.references',
                            colors=colors,
                            width=width,
                            height=height,
                            outfile=outfile,
                            export=export )
    }

    return( ranking )
}


#
# Subpathway level measures
#

.subLevelMeasures <- function( subpathway, type, topTerms=20 )
{
    supportedTerms <- .getExternalMeasures()

    if ( missing(type) )
    {
        type <- supportedTerms
    }

    failType <- type[!type %in% supportedTerms]
    if ( length(failType) > 0 )
    {
        message('Option(s) ', failType, ' not supported.')
    }

    # Access local Rdata diles
    targets <- .loadTermData( type )
    if (is.null(targets)) { return(NULL) }

    # Calculate p-values for each class
    N        <- nrow(targets)
    pValues  <- vector( mode='numeric', length=N )
    uTargets <- unique( as.vector( targets ) )
    res <- matrix(, nrow=N, ncol=5)
    for ( i in 1:N )
    {
        termTargetsAll   <- which(targets[i,] != '0')
        termTargetsOfSub <- which(!is.na(match(targets[i,], 
                            subpathway)))
        # Hypergeometric test compares sample to background
        # Success.sample, Success.bgd, Failure.bgd, Sample.size
        A <- length( termTargetsOfSub )
        B <- length( uTargets )
        C <- length( termTargetsAll )
        D <- length( sub )
        pValues[i] <- 1 - phyper(A, C, B-C, D)
        pValues[i] <- floor(pValues[i] * 10000)/1000
    }
    names(pValues) <- rownames(targets)


    # Keep significant terms
    idx <- which(pValues < 0.05)
    if ( length(idx) == 0 )
    {
        return(NULL)
    }
    pValues <- pValues[idx]
    targets <- targets[idx, , drop=FALSE]
    # Keep top terms
    idx <- order(pValues)
    pValues <- pValues[idx]
    targets <- targets[idx, , drop=FALSE]
    # Choose most significant classes
    if ( length(targets) == 0 ) { return(NULL) }  
    if ( length(pValues) < topTerms ) 
        { topTerms <- length(pValues) }
    pValues <- pValues[1:topTerms]
    targets <- targets[1:topTerms, , drop=FALSE]

    # Reduce target matrix with subpathway genes by removing rows and 
    # and columns not containing any subpathway genes
    idx <- targets %in% subpathway
    idx <- matrix( idx, nrow=nrow(targets) )*1
    targets[!idx] <- '0'

    # Remove rows with no subpathway genes
    idx <- targets == '0'
    if ( length(idx) > 0 )
    {
        idx <- rowSums(idx) != ncol(targets)
        targets <- targets[idx, , drop=FALSE]
    }
    # Remove columns with no subpathway genes
    idx <- targets == '0'
    if ( length(idx) > 0 )
    {
        idx <- colSums(idx) != nrow(targets)
        targets <- targets[, idx, drop=FALSE]
    }
    # Check if any data remains
    if ( length(targets) == 0 )
        { return(NULL) }  

    # Create edgelist for topmost results (class, subpathway genes)
    edgeList <- matrix(, nrow=length(targets), ncol=3)
    for ( i in seq_len(nrow(targets)) )
    {
        start <- 1 + (i-1)*ncol(targets)
        end   <- i*ncol(targets)
        genesPerTerm <- expand.grid(rownames(targets)[i], targets[i, ] )
        edgeList[start:end, 1:2] <- as.matrix(genesPerTerm)
        edgeList[start:end, 3] <- pValues[rownames(targets)[i]]
    }
    idx <- which(edgeList[, 2] != '0')
    edgeList <- edgeList[idx, , drop=FALSE]

    return ( edgeList )
}


subpathwayVisualization <- function( DEsubs.out, 
                                references='', 
                                submethod, 
                                subname,
                                colors=c('#FF0000FF', '#FF9900FF','#CCFF00FF',
                                        '#33FF00FF', '#00FF66FF', '#0066FFFF',
                                        '#3300FFFF', '#CC00FFFF', '#FF0099FF',
                                        '#EE82EEFF'), 
                                scale=1, 
                                shuffleColors=FALSE, 
                                outfiles='', 
                                export='pdf', 
                                verbose=TRUE )
{
    if ( 'all' %in% references )
    {
        references <- .getExternalMeasures()
        references <- c(references, 'GO', 'Disease', 'Regulator')
    }
    if ( length(scale) == 1 )
    {
        scale <- rep(scale, length(references))
    }
    if ( length(references) != length(scale) )
    {
        message('Number of references has to be equal to different 
            scaling optios. Setting scale to 1. ')
        scale <- rep(1, length(references))
    }
    if ( (!'' %in% outfiles)  && (length(outfiles) != length(references) ) )
    {
        message('Number of outfiles should be equal to number of references.')
        return(NULL)
    }
    if ( '' %in% outfiles ) { outfiles <- rep('', length(references) ) }

    out <- vector(mode='list', length=length(references))
    for ( i in seq_len(length(references)) )
    {
        out[[i]] <- .subpathwayVisualization(DEsubs.out=DEsubs.out,
                                reference=references[i],
                                submethod=submethod,
                                subname=subname,
                                colors=colors,
                                scale=scale[i],
                                shuffleColors=shuffleColors,
                                outfile=outfiles[i],
                                export=export,
                                verbose=verbose
                                )
    }
    names(out) <- references

    return( out )
}

.subpathwayVisualization <- function( DEsubs.out, reference, submethod, 
                                    subname, colors, scale, shuffleColors, 
                                    outfile, export, verbose)
{

    org <- DEsubs.out[['org']]
    edgeList <- DEsubs.out[['edgeList']]

    if ( org != 'hsa' ) 
        { return(NULL) } # Homo sapiens only
    if ( is.null(submethod) ) 
        { return(NULL) } 


    supportedMethods <- subpathwayTypes()
    supportedReferences <- .getExternalMeasures()
    supportedReferences <- c(supportedReferences, 'GO', 'Disease', 'Regulator')

    if (verbose)
    {
        message('Generating subpathway-level view ( ', reference,' )...', 
                                                        appendLF = FALSE)
    }

    unsupportedOptions <- submethod[!submethod %in% supportedMethods]
    if ( length(unsupportedOptions) > 0 )
    {
        message('Option(s) ', unsupportedOptions, ' not supported.')
        return(NULL)
    }
    unsupportedReferences <- reference[!reference %in% supportedReferences]
    if ( length(unsupportedReferences) > 0 )
    {
        message('Reference(s) ', unsupportedReferences, ' not supported.')
        return(NULL)
    }


    # Subpathway information extraction
    submethodName <- paste0('subAnalysis.', submethod)
    subpathway.nodelist <- DEsubs.out[[submethodName]][[subname]]

    # Subpathway (nodelist) to edgelist
    idx1 <- which(edgeList[, 1] %in% subpathway.nodelist)
    idx2 <- which(edgeList[, 2] %in% subpathway.nodelist)
    idx <- intersect(idx1, idx2)
    subpathway.edgelist <- edgeList[idx, ]

    # Change annotation from entrez to HGNC
    subpathway.nodelist <- .changeAnnotation(annData=subpathway.nodelist, 
                        org='hsa', choice='entrezToHGNC')
    subpathway.nodelist <- unname(subpathway.nodelist)


    # Special care for double and triple references
    references <- reference
    top <- 30
    agap <- 0; bgap <- 0; cgap <- 50
    if ( reference %in% 'GO' )
    {
        references <- c('GO_bp', 
                        'GO_cc', 
                        'GO_mf')
        agap <- 20; bgap <- 20; cgap <- 50
    }
    if ( reference %in% 'Disease' )
    {
        references <- c('Disease_OMIM', 'Disease_GAD')
        agap <- 20; bgap <- 0; cgap <- 50
    }
    if ( reference %in% 'Regulator' )
    {
        references <- c('miRNA', 'TF')
        agap <- 20; bgap <- 0; cgap <- 50
    }
    top <-  floor(top/length(references))

    edgeLists <- list()
    for ( i in seq_len(length(references)) )
    {
        subpathway.edgelist <- .subLevelMeasures(
                            subpathway=subpathway.nodelist, 
                            type=references[i], 
                            topTerms=top)
        if ( !is.null(subpathway.edgelist) )
        {
            if ( top > nrow(subpathway.edgelist) )
                { top <- nrow(subpathway.edgelist) }
            subpathway.edgelist <- subpathway.edgelist[1:top, , drop=FALSE]

            R <- list(subpathway.edgelist)
            names(R) <- references[i]
            edgeLists <- c(edgeLists,  R)
        }
    }

    if ( length(edgeLists) == 0 ) 
    { 
        if (verbose) { message('done.') }
        return(NULL) 
    }

    # Merge rows and columns indepedently for each edgelist.
    rowTerms <- c(); colTerms <- c()
    for ( i in seq_len(length(edgeLists)) )
    {
        if ( is.null(edgeLists[[i]]) ) { next() }
        rowTerms <- unique(c(rowTerms, edgeLists[[i]][, 1]))
        colTerms <- unique(c(colTerms, edgeLists[[i]][, 2]))
    }


    edgeList <- do.call(rbind, edgeLists)
    adjmat   <- matrix(0, nrow=length(rowTerms), ncol=length(colTerms))
    rownames(adjmat) <- rowTerms
    colnames(adjmat) <- colTerms
    adjmat.pval <- matrix(NA, nrow=length(rowTerms), ncol=length(colTerms))
    rownames(adjmat.pval) <- rowTerms
    colnames(adjmat.pval) <- colTerms    
    for ( i in seq_len(nrow(edgeList)) )
    {
        # For visualization purposes
        adjmat[ edgeList[i,1], edgeList[i,2] ] <- 1
        # For output
        adjmat.pval[edgeList[i,1], edgeList[i,2]] <- as.numeric(edgeList[i,3])
    }

    if ( shuffleColors  )
        { colors <- sample(colors, nrow(adjmat), replace=TRUE) }
    if ( !shuffleColors  )
        { colors <- rep(colors, length.out=nrow(adjmat)) }


    # The user has supplied a custom directory
    if ( (outfile != '') && ('pdf' %in% export) ) 
    { 
        .dir.create.rec(outfile) 
    }
    # No custom directory has been supplied
    if ( (outfile == '') && ('pdf' %in% export) )
    {
        dir.create(cache[['outDir']], showWarnings=FALSE, recursive=TRUE)
        outfile <- paste0(cache[['outDir']] , '//circos_', reference ,'.pdf')
    }

    .doCirclize( mat=adjmat, 
                colors=colors,
                agap=agap,
                bgap=bgap,
                cgap=cgap,
                degree=250,
                a=round(10/scale)/10,
                b=round(10/scale)/10,
                outfile=outfile,
                export=export)

    if (verbose) { message('done.') }

    return( adjmat.pval )
}


#
# Organism level measures
#
organismVisualization <- function( DEsubs.out, 
                            references='', 
                            topSubs=10, 
                            topTerms=20, 
                            colors=colorRampPalette(c("white", "red"))(100), 
                            export='pdf', 
                            width=7, 
                            height=6, 
                            outfiles='', 
                            verbose=TRUE )
{
    # Homo sapiens only
    org <- DEsubs.out[['org']]
    if (org != 'hsa') { return(NULL) }

    if ( 'all' %in% references )
    {
        references <- .getExternalMeasures()
    }

    if ( (!'' %in% outfiles)  && (length(outfiles) != length(references) ) )
    {
        message('Number of outfiles should be equal to number of references.')
        return(NULL)
    }
    if ( '' %in% outfiles ) { outfiles <- rep('', length(references) ) }

    out <- vector(mode='list', length=length(references))
    for ( i in seq_len(length(references)) )
    {
        out[[i]] <- .organismVisualization(DEsubs.out=DEsubs.out,
                                    references=references[i],
                                    topSubs=topSubs,
                                    topTerms=topTerms,
                                    colors=colors,
                                    export=export,
                                    width=width,
                                    height=height,
                                    outfile=outfiles[i],
                                    verbose=verbose)
    }
    names(out) <- references

    return( out )
}

.organismVisualization <- function( DEsubs.out, references, topSubs, topTerms,
                                    colors, export, width, height, 
                                    outfile, verbose )
{

    org <- DEsubs.out[['org']]
    type <- references
    supportedTerms <- .getExternalMeasures()

    if (verbose)
    {
        message(paste0('Generating organism-level view ( ', type , ' )...'), 
                appendLF = FALSE)
    }

    failType <- type[!type %in% supportedTerms]
    if ( length(failType) > 0 )
    { 
        message('Option(s) ', failType, ' not supported.')
        return(invisible())
    }

    subpathways <- DEsubs.out[[5]]
    termsPerSub <- list()
    pValuesPerTermPerSub <- list()
    pathNames <- c()


    top <- c(topSubs, topTerms)
    names(top) <- c('subpathways', 'terms') 
    
    if ( length(subpathways) < top['subpathways'] ) 
        { top['subpathways'] <- length(subpathways) }

    subpathways <- subpathways[1:top['subpathways']]
    subpathways <- subpathways[!is.na(subpathways)]
    offset <- 0.001

    # Find most significant terms for each subpathway
    for ( i in seq_len(length(subpathways)) )
    {
        sub <- unname(.changeAnnotation(annData=subpathways[[i]], org=org,
                    choice='entrezToHGNC'))
        edgeList <- .subLevelMeasures(subpathway=sub, type=type, 
                    topTerms=top['subpathways'])
        if ( !is.null(edgeList) )
        {
            termsPerSub[[i]] <- matrix( edgeList[, 1], nrow=1 )
            pval <- matrix( as.numeric(edgeList[, 3]) + offset, nrow=1 )
            pValuesPerTermPerSub[[i]] <- pval
            pathNames <- c(pathNames, names(subpathways)[i] )
        }
    }

    if ( length(termsPerSub) == 0) 
    {
        message('done', appendLF = TRUE) 
        return(NULL) 
    }
    
    # Convert lists to adjacency matrices
    names(termsPerSub) <- pathNames
    termsPerSub <- .unlistToMatrix(.fillMatrixList( termsPerSub ) )
    names(pValuesPerTermPerSub) <- pathNames
    pValuesPerTermPerSub <- .fillMatrixList( pValuesPerTermPerSub )
    pValuesPerTermPerSub <- .unlistToMatrix( pValuesPerTermPerSub )

    # Find topmost terms amongst the terms for all subpathways
    uniqueTerms <- unique(as.vector(termsPerSub))
    uniqueTerms <- uniqueTerms[uniqueTerms != '0']
    counts <- vector(mode='numeric', length=length(uniqueTerms))
    for ( i in seq_len(length(uniqueTerms)) )
    {
        counts[i] <- length(which(as.vector(termsPerSub) == uniqueTerms[i]))
    }
    names(counts) <- uniqueTerms
    topTerms <- names(head(sort(counts, decreasing=TRUE), top['terms']))

    # Keep only topmost terms for each subpathway
    idx <- termsPerSub %in% topTerms
    if ( length(idx) > 0 ) { termsPerSub[!idx] <- NA }

    # Create an edgelist of topmost terms for all subpathways
    termsPerSub.edgeList <- matrix(, nrow=length(termsPerSub), ncol=3)
    for ( i in seq_len(nrow(termsPerSub)) )
    {
        if ( is.na(rownames(termsPerSub)[i]) ) { next() }

        start <- 1 + (i-1)*ncol(termsPerSub)
        end   <- i*ncol(termsPerSub)
        data1 <- expand.grid( rownames(termsPerSub)[i], termsPerSub[i, ] )
        data2 <- expand.grid(   rownames(termsPerSub)[i], 
                                pValuesPerTermPerSub[i, ] - offset )

        data.all <- cbind(data1, data2[, 2])
        termsPerSub.edgeList[start:end, ] <- as.matrix(data.all)        
    }

    # Remove NA's
    idx <- which(is.na(termsPerSub.edgeList[, 2]))
    if ( length(idx) > 0 ) 
    { 
        termsPerSub.edgeList <- termsPerSub.edgeList[-idx, , drop=FALSE] 
    }

    termsPerSub.df <- as.data.frame(termsPerSub.edgeList, 
                                                    stringsAsFactors=FALSE)
    colnames(termsPerSub.df) <-  c('Subpathway', 'Term', 'pValue')

    termsPerSub.df[, 'pValue'] <- as.numeric(termsPerSub.df[, 'pValue'])


    .matrixVisualization(dataVis=termsPerSub.df, type='dotplot', 
                        title=c('Subpathway', type, 'Q-values'),
                        colors=colors,
                        width=width, height=height,
                        outfile=outfile,
                        export=export)

    if (verbose)
        { message('done', appendLF = TRUE) }
    

    return( termsPerSub.df )
}


