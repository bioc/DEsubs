.fillMatrixList  <- function( limat, maxLen )
{
    # Zero pads a list of matrices (limat) each having different number of 
    # columns. If the maximum number of columns is availiable, providing 
    # is as an argument speeds up computation time (maxlen).
    if (is.null(limat))     { return(limat) }

    # If input is not a list, return original data
    if (!is.list(limat))    { return(limat) }

    lens <- c()
    for (i in seq_len(length(limat)) )
    {    
        if (!is.null(limat[[i]]))
        {
            lens <- c(lens, .getLengths(limat[[i]]))
        }
    }

    if (is.null(lens)) { return(NULL) }


    # Set number of columns.
    if (missing(maxLen)) { maxLen <- max(lens, na.rm=TRUE) }


    # If no filling is necessary, return original data
    if (length(unique(lens)) == 1) { return(limat) }

    # Enforce maximum length on each subpath of each pathway
    res <- limat
    for (i in seq_len(length(limat)) )
    {
        submat <- limat[[i]]
        if (is.null(submat)) { next() }
        replmat <- matrix(0, nrow=nrow(submat), ncol=maxLen)
        for (j in seq_len(nrow(submat)) )
        {   
            replmat[j,] <- c(submat[j,], rep(0, maxLen - length(submat[j,])))
        }
        res[[i]] <- replmat
    }

    return(res)
}

.unlistToMatrix  <- function( limat, mode='rbind' )
{
    # Reshapes a list of matrices to one matrix in a full vectorized fashion
    if (is.null(limat))     { return(limat) }

    # If input is already a matrix, return original data.
    if (class(limat) == 'matrix') { return(limat) }

    # Find what type of data the list holds
    type <- NULL
    for (i in seq_len(length(limat)) )
    {
        if(!is.null(limat[[i]]))
        {
            if (is.vector(limat[[i]])) type <- 'vector'
            if (is.matrix(limat[[i]])) type <- 'matrix'
            break()
        }
    }

    if (is.null(type)) { return(NULL) }

    # List of vectors of variable size
    if (type == 'vector')
    {
        lens <- sapply(limat, function(x) { length(x) })
        M    <- matrix(0, nrow=length(limat), ncol=max(lens))
        for (i in seq_len(length(limat)) ) 
        {
            M[i,1:lens[i]] <- limat[[i]] 
        }
        rownames(M) <- names(limat)
    }

    # A list of matrices with at least one common dimension size
    # rbind:same number of columns, cbind:same number of rows
    if (type == 'matrix')
    {
        M      <- do.call(mode, limat)
        lnames <- names(limat) 

        if (length(lnames) > 0)
        {
            lens   <- sapply(limat, function(x) 
                                { 
                                    if (!is.null(x)) nrow(x) else 0 
                                } )
            rnames <- vector(mode='numeric', length=nrow(M))
            ctr   <- 0
            for (i in seq_len(length(lnames)) )
            {
                if (lens[i] > 0)
                {
                    idx         <- (ctr+1) : (ctr<-ctr+lens[i])
                    rnames[idx] <- rep(lnames[i], lens[i])               
                }
            }
            rownames(M) <- rnames
        }
    }
    return(M)
}

.getLengths      <- function( mat )
{
    # 
    # Count length of non zero elements in each row of a matrix.
    # Each row must consist of a prefix with non zero elements, 
    # and a suffix of consecutive zeros denote absence of data.
    # Ideal for matrices with a large number of rows.
    #
    esub <- cbind(mat, rep(0, nrow(mat)))
    esub <- rbind(esub, c(-1, rep(0, ncol(mat))))
    df   <- which(diff(which(c(t(esub)) != 0)) - 1 > 0)
    len  <- c(df[1], diff(df))

    return(len)
}


.dir.create.rec <- function(file) 
{
    if( !file.exists(file) ) 
    {
        .dir.create.rec(dirname(file))
        if ( gsub('\\.', '', file) == file  ) 
        {
            dir.create(file, recursive=FALSE, showWarnings=TRUE) 
        }
    }
} 

