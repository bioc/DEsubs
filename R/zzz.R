cache <- new.env()


.onLoad <- function(libname, pkgname)
{

    #
}


.onAttach <- function(libname, pkgname)
{
    # Set the default location of the base directory
    path <- switch(.Platform[['OS.type']], unix = path.expand("~"),
                    windows= file.path(gsub("\\\\", "/",
                    Sys.getenv("USERPROFILE")), "AppData"))
    opt <- getOption("DEsubs_CACHE", paste0(path, '/DEsubs'))
    baseDir <- Sys.getenv("DEsubs_CACHE", opt)

    cache[['baseDir']] <- baseDir
    cache[['outDir']] <- paste0(baseDir, '//Output')
    cache[['datDir']] <- paste0(baseDir, '//Data')
}


.onUnload <- function(libname) 
{

    # 
}

