
subpathwayTypes <- function(grouping='all')
{

    subTypes <- .subpathwayTypes()

    supportedMethods <- c(  subTypes[['subStreamCases']], 
                            subTypes[['neighborCases']], 
                            subTypes[['linearCases']], 
                            subTypes[['communityCases']], 
                            subTypes[['componentCases']])

    if ( is.null(grouping) || grouping == '' ) { return(NULL) }


    if ( length(grouping) == 1 && grepl('^all', grouping) )
    {
        subType <- supportedMethods
        type <- gsub('all.', '', grouping)
        intMeasures <- c('bwd', 'fwd', 'stream', 'neighbourhood', 
                        'cascade', 'community', 'component', 
                        'topological', 'functional', 'DEG')
        extMeasures <- .getExternalMeasures()

        if ( type %in% c(intMeasures, extMeasures) )
        {
            supportedMethods <- subType[grepl(type, subType)]
        }
    }
    if ( length(grouping) == 1 && !grepl('^all', grouping) )
    {
        supportedMethods <- grouping
    }

    return( supportedMethods )
}


.subpathwayTypes <- function()
{
    stream.choices <- c('fwd',
                        'bwd')
    commun.choices <- c('walktrap', 
                        'edge_betweenness', 
                        'fast_greedy', 
                        'leading_eigen', 
                        'infomap', 
                        'louvain'
                        )
    compon.choices <- c('decompose', 
                        'max_cliques',
                        'cliques',
                        'coreness'
                        )
    topological.choices <- c(
                        'topological.degree', 
                        'topological.betweenness', 
                        'topological.closeness', 
                        'topological.hub_score', 
                        'topological.eccentricity', 
                        'topological.page_rank',
                        'topological.start_nodes')

    functional.choices <- paste0('functional.', .getFunctionalMeasures())

    source.choices <- c(topological.choices, functional.choices)

    cases <- expand.grid( stream.choices, 'neighbourhood', source.choices)
    neighborCases <- apply(cases, 1, function(x) {paste0(x, collapse='.')})

    cases <- expand.grid( stream.choices, 'stream', source.choices)
    subStreamCases <- apply(cases, 1, function(x) {paste0(x, collapse='.')})

    cases <- expand.grid( stream.choices, 'cascade', source.choices)
    linearCases <- apply(cases, 1, function(x) {paste0(x, collapse='.')})

    cases <- expand.grid( 'community', commun.choices)
    communityCases <- apply(cases, 1, function(x) {paste0(x, collapse='.')})

    cases <- expand.grid( 'component', compon.choices[1:2])
    compCases_a <- apply(cases, 1, function(x) {paste0(x, collapse='.')})
    cases <- expand.grid( 'component.', paste0(3:9, '-', compon.choices[3]) )
    compCases_b <- apply(cases, 1, function(x) {paste0(x, collapse='')})
    cases <- expand.grid( 'component.', paste0(3:9, '-', compon.choices[4]) )
    compCases_c <- apply(cases, 1, function(x) {paste0(x, collapse='')})
    componentCases <- c(compCases_a, compCases_b, compCases_c)

    return(list('neighborCases'=neighborCases,
                'subStreamCases'=subStreamCases,
                'linearCases'=linearCases,
                'communityCases'=communityCases,
                'componentCases'=componentCases))
}


.subpathwayAnalysis <- function( edgeList, method=c(), DEgenes, a=3, b=10, 
                                org, verbose=TRUE )
{
    # a: Number of neighbors
    # b: Number of top-nodes 

    supportedMethods <- subpathwayTypes()

    unsupportedOptions <- method[!method %in% supportedMethods]
    if ( length(unsupportedOptions) > 0 )
    {
        message('Option(s) ', unsupportedOptions, ' not supported.')
        out        <- list(NULL)
        names(out) <- paste0('subAnalysis.', method)
        return( out )
    }

    if (verbose)
    {
        message('Performing subpathway analysis ( ', method , ' )...', 
                appendLF = FALSE) 
    }

    # Create an igraph object from an edgelist
    gi <- graph_from_edgelist(as.matrix(edgeList[, 1:2]))

    # Change direction of substream if necessary
    if ( gsub('bwd.', '', method) != method  )
    {
        edgeList <- get.edgelist(gi)
        gi <- graph_from_edgelist(cbind(edgeList[, 2], edgeList[, 1]))
    }
    
    if ( nrow(edgeList) == 0 )
    {
        if (verbose) { message('done.') }
        return (NULL)
    }

    subTypes <- .subpathwayTypes()

    if ( method %in% subTypes[['neighborCases']] )
    {
        sourceMeasure <- gsub('fwd.neighbourhood.', '', method)
        sourceMeasure <- gsub('bwd.neighbourhood.', '', sourceMeasure)

        nodes <- .measureToNodes(graph=gi,
                                measure=sourceMeasure,
                                org=org,
                                DEgenes=DEgenes )
        sourceNodes <- nodes[1 : min(b, length(nodes))]
        
        if (length(nodes) == 0)
        { 
            out        <- list(NULL)
            names(out) <- paste0('subAnalysis.', method)
            if (verbose) { message('done.') }
            return(out)
        }

        R <- vector(mode='list', length=length(sourceNodes))
        for (i in seq_len(length(sourceNodes)) )
        {
            R[[i]] <- names(unlist(ego(graph=gi, order=a, 
                                    nodes=sourceNodes[i])))
        }
        names(R) <- paste0('sub', seq_len(length(R)))

        out        <- list(R)
        names(out) <- paste0('subAnalysis.', method)
    }
    if ( method %in% subTypes[['subStreamCases']] )
    {
        sourceMeasure <- gsub('fwd.stream.', '', method)
        sourceMeasure <- gsub('bwd.stream.', '', sourceMeasure)
        nodes <- .measureToNodes(graph=gi,
                                measure=sourceMeasure,
                                org=org,
                                DEgenes=DEgenes )
        sourceNodes <- nodes[1 : min(b, length(nodes))]
        
        if (length(nodes) == 0)
        { 
            out        <- list(NULL)
            names(out) <- paste0('subAnalysis.', method)
            if (verbose) { message('done.') }
            return(out)
        }

        R <- vector(mode='list', length=length(sourceNodes))
        for (i in seq_len(length(sourceNodes)) )
        {
            res <- names(dfs(graph=gi, root=sourceNodes[i], 
                            unreachable=FALSE)[['order']])
            R[[i]] <- res[which(!is.na(res))]
        }
        names(R) <- sourceNodes

        # Keep subpathways with more than two members
        R <- R[which( lapply(R, function(x) { length(x) }) >= 3 )]
        if ( length(R) > 0 ) { names(R) <- paste0('sub', seq_len(length(R)) ) }

        out        <- list(R)
        names(out) <- paste0('subAnalysis.', method)
    }
    if ( method %in% subTypes[['linearCases']] )
    {
        sourceMeasure <- gsub('fwd.cascade.', '', method)
        sourceMeasure <- gsub('bwd.cascade.', '', sourceMeasure)

        # Find all source nodes
        nodes <- .measureToNodes(graph=gi,
                                measure=sourceMeasure,
                                org=org,
                                DEgenes=DEgenes )
        sourceNodes <- nodes[1 : min(b, length(nodes))]

        if (length(nodes) == 0)
        { 
            out        <- list(NULL)
            names(out) <- paste0('subAnalysis.', method)
            if (verbose) { message('done.') }
            return(out)
        }

        # Simplify graph
        adjmat <- get.adjacency(gi)
        gi <- graph.adjacency(triu(adjmat))

        # Find all destination nodes
        destinNodes <- names(sort(which(degree(gi, mode='out', 
                                loops=FALSE) == 0), decreasing=TRUE))
        sourceNodes <- sourceNodes[!sourceNodes %in% destinNodes]
        destinNodes <- destinNodes[!destinNodes %in% sourceNodes]
                
        # Find all linear subpathways between each pair of start/end nodes
        N <- length(sourceNodes)*length(destinNodes)
        lpaths <- vector(mode='list', length=N)
        lnames <- vector(mode='numeric', length=N)
        cnt <- 1

        lpaths <- NULL
        if ( N > 0 )
        {
            for ( sourceNode in sourceNodes )
            {
                for ( destinNode in destinNodes )
                {
                    lpaths[[cnt]] <- all_simple_paths( gi, 
                                                from = sourceNode, 
                                                to = destinNode,
                                                mode="out")
                    lnames[cnt] <- paste(sourceNode, destinNode, sep='-')
                    cnt <- cnt + 1
                }
            }

            combinations <- expand.grid( sourceNodes, destinNodes )
            names(lpaths) <- lnames
            names(lpaths) <- paste0('sub', seq_len(length(lpaths)) )

            # Remove NULL results and flatten first level of results 
            lpaths <- do.call(c, lpaths)
            # Keep subpathways with more than two members
            lpaths <- lpaths[which(lapply(lpaths, length) > 2)]
            # Keep genes ids
            lpaths <- lapply(lpaths, function(x) { names(x) } )
        }

        out        <- list(lpaths)
        names(out) <- paste0('subAnalysis.', method)
    }
    if ( method %in% subTypes[['communityCases']] )
    {
        communityHandlers <- c( 
            'walktrap'=function(x) 
                            { cluster_walktrap(x) },
            'edge_betweenness'=function(x) 
                            { cluster_edge_betweenness(x) },
            'fast_greedy'=function(x) 
                            { cluster_fast_greedy(as.undirected(x)) },
            'leading_eigen'=function(x) 
                            { cluster_leading_eigen(as.undirected(x)) },
            'infomap'=function(x) 
                            { cluster_infomap(as.undirected(x)) },
            'louvain'=function(x) 
                            { cluster_louvain(as.undirected(x)) })

        gr <- communityHandlers[[gsub('community.', '', method)]](gi)
        R <- vector( mode='list', length=length(gr) )
        if ( length(gr) > 0 )
            { for ( i in seq_len(length(gr)) ) { R[[i]] <- (gr[[i]]) } }  
        if ( length(R) > 0 ) 
            { names(R) <- paste0('sub', seq_len(length(R)) ) }

        out <- list(R)
        names(out) <- paste0('subAnalysis.', method)
    }
    if ( method %in% subTypes[['componentCases']] )
    {
        if ( method == 'component.max_cliques' )
        { 
            gr <- max_cliques(as.undirected(gi), a, b)
            R <- vector( mode='list', length=length(gr) )
            if ( length(gr) > 0 )
            {
                for ( i in seq_len(length(gr)) ) 
                                { R[[i]] <- names(unlist(gr[i])) }
            }
        }

        if ( method == 'component.decompose')
        {
            gr  <- decompose(gi)
            R   <- lapply(gr, function(x) { as_data_frame(x, what="edges") } )
        }

        if ( gsub('-cliques', '', method) !=  method )
        {
            p <- gsub('-cliques', '', method)
            k <- as.numeric(gsub('component.', '', p))

            g <- as_graphnel( gi )
            ksubs <- kCliques(ugraph(g))[paste0(k, '-cliques')][[1]]
            ksubs <- lapply(ksubs, function(x) { matrix(x, nrow=1) } )

            # Extract valid pairs
            edgeList <- as_data_frame(gi, what="edges")
            
            # Consider only actual interactions between genes.
            if ( length(ksubs) > 0 )
            {
                ksubs <- .unlistToMatrix(.fillMatrixList(ksubs))
                R <- vector(mode='list', length=nrow(ksubs))
                for (j in seq_len(nrow(ksubs)) )
                {
                    # Edge list indexes
                    r1  <- as.numeric(is.element( edgeList[,1], ksubs[j, ]) )
                    r2  <- as.numeric(is.element( edgeList[,2], ksubs[j, ]) )
                    idx <- which((r1 + r2) == 2)
                    if ( length(idx) > 0 )
                    { 
                        R[[j]] <- unique(as.vector(t(edgeList[idx, ])))
                    }
                }
            }else{ R <- NULL}

        }
        if ( gsub('-coreness', '', method) !=  method )
        {
            p <- gsub('-coreness', '', method)
            p <- as.numeric(gsub('component.', '', p))

            corDist  <- coreness(gi)
            R <- names(corDist)[which(corDist == p)]
        }
        if ( length(R) > 0 ) { names(R) <- paste0('sub', seq_len(length(R)) ) }


        out <- list(R)
        names(out) <- paste0('subAnalysis.', method)
    }
    
    if (verbose) { message('done.') }


    return(out)
}


.getFunctionalNodes <- function( graph, targets, org )
{
    # Extract nodes from graph and keep term related ones
    graphGenes <- names(V(graph))
    # Unique target genes in entrez 
    uGenesFromTargets <- unique(as.vector(targets))
    uGenesFromTargets <- unname(.changeAnnotation(annData=uGenesFromTargets, 
                                org='hsa', choice='HGNCtoEntrez'))
    uGenesFromTargets <- uGenesFromTargets[!is.na(uGenesFromTargets)]
    # Keep targets intersecting with graph genes
    targetGenesInGraph <- graphGenes[graphGenes %in% uGenesFromTargets]
    targetGenesInGraph <- unname(.changeAnnotation(annData=targetGenesInGraph, 
                                org='hsa', choice='entrezToHGNC'))

    # Filter target genes with graph genes
    idx <- matrix(targets %in% targetGenesInGraph, nrow=nrow(targets))*1
    targets[!idx] <- NA

    nodes <- names(sort(table(targets), decreasing=TRUE))
    # Change to entrez gene annotation
    nodes <- unname(.changeAnnotation(annData=nodes, org=org, 
                    choice='HGNCtoEntrez'))

    return (nodes)
}


.measureToNodes <- function ( graph, measure, org, DEgenes=NULL )
{
    topologicalHandlers <- .getTopologicalHandlers()

    if ( grepl('topological', measure) )
    {
        func <- topologicalHandlers[[gsub('topological.', '', measure)]]
        nodes <- names(sort(func(graph), decreasing=TRUE)) 
    }

    if ( measure == 'functional.DEG' ) 
        { nodes <- names(sort(DEgenes, decreasing=FALSE)) }

    if ( measure %in% paste0('functional.', .getExternalMeasures() ) )
    {
        if ( org != 'hsa' )
        { 
            out <- list(NULL)
            names(out) <- paste0('subAnalysis.', measure)
            return( out )
        }

        # Find genes with most occurences in the selected term
        measure <- gsub('functional.', '', measure)
        targets <- .loadTermData( type=measure )
        if ( is.null(targets) )
            { return(NULL) }
        nodes <- .getFunctionalNodes(graph=graph, targets=targets, org=org)
    }

    return( nodes )
}

.getTopologicalHandlers <- function()
{
    topologicalHandlers <- c( 
                'degree'=function(x) { degree(x) },
                'betweenness'=function(x) { betweenness(x) },
                'closeness'=function(x) { closeness(x) },
                'hub_score'=function(x) { hub_score(x)[['vector']] },
                'eccentricity'=function(x) { eccentricity(x) },
                'page_rank'=function(x) { page_rank(x)[['vector']] },
                'start_nodes'=function(x) 
                    { which(degree(x, mode='in', loops=FALSE) == 0) } )

    return( topologicalHandlers )
}

#
# Data-related functions
#
.changeAnnotation <- function(annData, org, choice)
{
    if ( org == 'hsa' )
    {
        dir <- system.file('extdata//Data', package='DEsubs')

        load(paste(dir, 'libraryEntrezToHGNC.RData', sep='//'), 
            e <- new.env())
        libraryEntrezToHGNC <- e[['libraryEntrezToHGNC']]
        if ( choice == 'entrezToHGNC' )
        {
            annData <- libraryEntrezToHGNC[annData]    
        }
        if ( choice == 'HGNCtoEntrez' )
        {
            libraryHGNCtoEntrez <- names(libraryEntrezToHGNC)
            names(libraryHGNCtoEntrez) <- libraryEntrezToHGNC
            annData <- libraryHGNCtoEntrez[annData]    
        }
    }

    return(annData)
}

.getExternalMeasures <- function()
{
    defaultReferences <- .getDefaultReferences()

    # The default datasets for external measures are stored in the Data
    # folder within the package directory. The user however can include 
    # any number of gene-sets within cache[['datDir']] directory.
    # The availiable external measures will be the union of gene-sets in
    # both directories. 

    files <- list.files(system.file('extdata//Data', package='DEsubs'))
    otherFiles <- c('libraryEntrezToHGNC.RData', 
                    'edgeLists.RData',
                    'libraryEntrezToExternalNomenclature.RData')
    files <- files[-which(files %in% otherFiles)]

    # Find if any of the default files have been deleted
    idx <- which(!paste0(defaultReferences, '.RData') %in% files)
    if ( length(idx) > 0 )
    {
        message( 'References ', paste0(files[idx], collapse=', '), 
                    ' missing.')
    }
    supportedReferences <- gsub('.RData', '', files)

    # Search for new gene sets
    userReferences <- gsub('.RData', '', list.files(cache[['datDir']]))

    supportedReferences <- unique(c(supportedReferences, userReferences))


    return( supportedReferences )
}

.getFunctionalMeasures <- function()
{

    return( c('DEG', .getExternalMeasures()) )
}

.getDefaultReferences <- function()
{
    # Default external references stored within the package library
    defaultReferences <- c( 'KEGG',
                            'GO_bp',
                            'GO_cc',
                            'GO_mf',
                            'Disease_OMIM',
                            'Disease_GAD',
                            'Drug_DrugBank',
                            'miRNA',
                            'TF')
    return( defaultReferences )
}

#
# Base data
#

.loadTermData <- function( type )
{
    references <- .getExternalMeasures()
    defaultReferences <- .getDefaultReferences()
    userReferences <- references[!references %in% defaultReferences ]

    if ( type %in% userReferences )
        { dir <- cache[['datDir']] }
    if ( type %in% defaultReferences )
        { dir <- system.file('extdata//Data', package='DEsubs') }

    if ( type %in% references )
    {
        iFile <- paste0( dir, '//' ,type, '.RData' )
        load(iFile, e <- new.env())        
        targetsPerClass <- e[['targetsPerClass']]
    } else{ message('Type ', type, ' not supported.') }

    if ( is.null(targetsPerClass) )
        { message('Set ', type, ' does not contain valid targers.') }

    return( targetsPerClass )
} 



