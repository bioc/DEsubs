DEsubs <- function( org, mRNAexpr, mRNAnomenclature, pathways, 
                    DEtool, DEpar, CORtool, CORpar, subpathwayType,
                    rankedList=NULL, verbose=TRUE)
{
    # Create the adjacency matrix
    dataNet <- .constructNetwork(org=org, 
                            mRNAexpr=mRNAexpr, 
                            mRNAnomenclature=mRNAnomenclature, 
                            pathways=pathways)
    mRNAexpr <- dataNet[['mRNAexpr']]
    edgeList <- dataNet[['edgeList']]

    # Filter the nodes and the edges of the organism's pathways network
    lens     <- suppressWarnings(split(1:ncol(mRNAexpr), 1:2))
    lens     <- unname(sapply(lens, function(x) { length(x) }))
    net      <- .pruneNetwork(  edgeList=edgeList,
                                mRNAexpr=mRNAexpr, 
                                DEGchoice=DEtool, 
                                DEGthresh=DEpar, 
                                classes=c(rep(1, lens[1]), rep(2, lens[2]) ),
                                corr_threshold=CORpar, org=org,
                                CORtool=CORtool,
                                rankedList=rankedList,
                                verbose=verbose)

    edgeList <- net[['edgeList']]
    DEgenes  <- net[['DEgenes']]
    output   <- list('org'=org, 'mRNAnomenclature'=mRNAnomenclature,
                    'edgeList'=edgeList, 'DEgenes'=DEgenes)

    # Choose specific subpathways
    subTypes <- subpathwayTypes(grouping=subpathwayType)

    # Subpathway analysis
    for (subType in subTypes)
    {
        subs <- .subpathwayAnalysis( edgeList=edgeList, 
                                    method=subType,
                                    DEgenes=DEgenes,
                                    org=org,
                                    verbose=verbose )
        output <- c(output, subs)
    }

    return( output )
}


.DEanalysis <- function( count.matrix, DEGchoice, classes )
{
    # 
    # Differential expression analysis using various DE analysis tools from
    # 
    # Soneson,C. and Delorenzi,M. (2013) A comparison of methods for 
    # differential expression analysis of RNA-seq data. BMC bioinformatics, 
    # 14(1), 1.
    #

    if ( missing(count.matrix) ) { message('Please supply a matrix.') }
    if ( missing(classes) )      { message('Please supply the classes.') }
    if ( missing(DEGchoice) )    { message('Please supply a type.') }

    supportedMethods <- c(
        'edgeR', 'DESeq2', 'EBSeq', 'NBPSeq', 'voom+limma', 'vst2+limma', 'TSPM')

    if ( DEGchoice == 'edgeR' )
    {
        # run edgeR
        edgeR.dgelist <- DGEList(counts=count.matrix,  group=factor(classes))
        edgeR.dgelist <- calcNormFactors(edgeR.dgelist,  method="TMM")
        edgeR.dgelist <- estimateCommonDisp(edgeR.dgelist)
        edgeR.dgelist <- estimateTagwiseDisp(edgeR.dgelist, trend="movingave")
        edgeR.test    <- exactTest(edgeR.dgelist)
        edgeR.pvalues <- edgeR.test[['table']][['PValue']]
        genes             <- rownames(edgeR.test[['table']])
        edgeR.adjpvalues  <- p.adjust(edgeR.pvalues,  method="BH")
        adjpvalues        <- edgeR.adjpvalues
        names(adjpvalues) <- genes

        return(adjpvalues)
    }
    if ( DEGchoice == 'voom+limma' )
    {
        # voom+limma
        nf                <- calcNormFactors(count.matrix,  method="TMM")
        voom.data         <- voom(  count.matrix, 
                                    design=model.matrix(~factor(classes)), 
                                    lib.size=colSums(count.matrix)*nf)
        voom.data[['genes']] <- rownames(count.matrix)
        voom.fitlimma        <- lmFit( voom.data, 
                                    design=model.matrix(~factor(classes)))
        voom.fitbayes     <- eBayes(voom.fitlimma)
        voom.pvalues      <- voom.fitbayes[['p.value']][, 2]
        voom.adjpvalues   <- p.adjust(voom.pvalues, method="BH")
        voom.genes        <- rownames(voom.fitbayes[['p.value']])
        adjpvalues        <- voom.adjpvalues
        names(adjpvalues) <- voom.genes

        return(adjpvalues)
    }
    if ( DEGchoice == 'EBSeq' )
    {
        # run EBSeq
        sizes       <-  MedianNorm(count.matrix)
        EBSeq.test  <-  suppressMessages( EBTest(Data=count.matrix, 
                        Conditions=factor(classes), sizeFactors=sizes,  
                        maxround=10))

        EBSeq.ppmat             <-  GetPPMat(EBSeq.test)
        EBSeq.probabilities.DE  <-  EBSeq.ppmat[, "PPDE"]
        EBSeq.lFDR  <-  1 - EBSeq.ppmat[, "PPDE"]
        EBSeq.FDR   <-  rep(NA,  length(EBSeq.lFDR))
        names(EBSeq.FDR) <- names(EBSeq.lFDR)
        for  (i  in  seq_len(length(EBSeq.lFDR)) )  
        {
            idx          <- which(EBSeq.lFDR <= EBSeq.lFDR[i])
            EBSeq.FDR[i] <- mean(EBSeq.lFDR[idx])
        }
        adjpvalues <- EBSeq.FDR

        return(adjpvalues)
    }
    if ( DEGchoice == 'NBPSeq' )
    {
        # NBPSeq
        NBPSeq.dgelist      <- DGEList( counts=count.matrix, 
                                        group=factor(classes))
        NBPSeq.dgelist      <- calcNormFactors(NBPSeq.dgelist, method="TMM")
        NBPSeq.norm.factors <- as.vector(
                            NBPSeq.dgelist[['samples']][['norm.factors']])
        
        capture.output(NBPSeq.test <- nbp.test(counts=count.matrix, 
            grp.ids=classes, grp1=1,grp2=2, norm.factors=NBPSeq.norm.factors))
        NBPSeq.pvalues    <- NBPSeq.test[['p.values']]
        NBPSeq.adjpvalues <- NBPSeq.test[['q.values']]
        adjpvalues        <- NBPSeq.adjpvalues 
        names(adjpvalues) <- rownames(NBPSeq.test[['counts']])

        return(adjpvalues)
    }
    if ( DEGchoice == 'TSPM' )
    {
        TSPM.dgelist <- DGEList(counts = count.matrix, group = factor(classes))
        TSPM.dgelist  <- calcNormFactors(TSPM.dgelist, method = "TMM")
        v1 <- as.vector(TSPM.dgelist[['samples']][['norm.factors']])
        v2 <- as.vector(TSPM.dgelist[['samples']][['lib.size']])
        norm.lib.sizes  <- v1 * v2 
        TSPM.test <- TSPM(counts = count.matrix, x1 = factor(classes), 
                    x0 = rep(1, length(classes)), lib.size = norm.lib.sizes)
        TSPM.pvalues    <- TSPM.test[['pvalues']]
        TSPM.adjpvalues <- TSPM.test[['padj']]
        adjpvalues      <- TSPM.adjpvalues
        genes           <- rownames(TSPM.dgelist[['counts']])
        names(adjpvalues) <- genes

        return(adjpvalues)
    }
    if ( DEGchoice == 'DESeq2' )
    {
        # run DESeq2
        sink(tempfile())
        DESeq2.ds <- DESeq2::DESeqDataSetFromMatrix(countData = count.matrix,
                        colData = data.frame(condition = factor(classes)), 
                        design = ~ condition)
        DESeq2.ds <- DESeq2::estimateSizeFactors( DESeq2.ds )
        DESeq2.ds <- DESeq2::estimateDispersions( DESeq2.ds, fitType="local")
        DESeq2.ds <- DESeq2::nbinomWaldTest( DESeq2.ds )
        DESeq2.results <- DESeq2::results(DESeq2.ds )
        adjpvalues <- p.adjust(DESeq2.results[['pvalue']], method="BH")
        names(adjpvalues) <- rownames(count.matrix)
        sink()

        return(adjpvalues)
    }
    if ( DEGchoice == 'vst2+limma' )
    {
        # vst(DESeq2)+limma
        
        sink(tempfile())
        DESeq2.ds <- DESeq2::DESeqDataSetFromMatrix(countData = count.matrix,
                        colData = data.frame(condition = factor(classes)), 
                        design = ~ condition)
        DESeq2.ds <- DESeq2::estimateSizeFactors( DESeq2.ds )
        DESeq2.ds <- DESeq2::estimateDispersions( DESeq2.ds, fitType="local")
        DESeq2.vst <- DESeq2::getVarianceStabilizedData( DESeq2.ds )

        DESeq2.vst.fitlimma <- lmFit(  DESeq2.vst,
                                        design=model.matrix(~factor(classes)))
        DESeq2.vst.fitbayes <- eBayes(DESeq2.vst.fitlimma)
        DESeq2.vst.pvalues  <- DESeq2.vst.fitbayes[['p.value']][, 2]
        genes               <- rownames(DESeq2.vst.fitbayes[['p.value']])
        DESeq2.vst.adjpvalues <- p.adjust(DESeq2.vst.pvalues, method="BH")
        adjpvalues          <- DESeq2.vst.adjpvalues 
        names(adjpvalues)   <- genes
        sink()

        return(adjpvalues)
    }

    if  ( !is.null(DEGchoice) )
    {
        unsupportedOptions <- DEGchoice[!DEGchoice %in% supportedMethods]
        if ( length(unsupportedOptions) > 0 )
        {
            message('Option ', unsupportedOptions, ' not supported.')
            message('Supported options are ', 
                        paste0(supportedMethods, collapse=', '), '.') 
            return( NULL ) 
        }
    }



    return( NULL )
}

