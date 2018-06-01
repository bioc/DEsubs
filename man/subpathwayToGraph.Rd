

\name{subpathwayToGraph}
\alias{subpathwayToGraph}
\title{Subpathway To Graph}
\usage{
subpathwayToGraph( DEsubs.out, submethod, subname, colors, size, export, 
                    width, height, outfile, verbose )
}
\arguments{
\item{DEsubs.out}{Return value from \code{\link{DEsubs}}}
\item{submethod}{Subpathway extraction type selection 
    (get all availiable options from \link{subpathwayTypes})}
\item{subname}{Subpathway name as contained in \code{\link{DEsubs}} return 
value, i.e. 'sub1'. }
\item{colors}{ A custom color mode which overrides the default settings. }
\item{size}{ A vector storing width and height of the barplot
visualization. }
\item{export}{ A set of options for exporting subpathway data. Possible options
are 'plot', 'pdf', 'edgelist', 'json', 'gml', 'ncol', 'lgl', 'graphml','dot'.}
\item{width}{ The width of the printable area (pdf) }
\item{height}{ The height of the printable area (pdf) }
\item{outfile}{ Output file name of the visualization. If multiple export
types have been selected, the outfile should have an extension '.*'}
\item{verbose}{ TRUE to display informative messages, FALSE to hide. }
}
\value{
No value is returned.
}
\description{
Subpathway plotting as a graph.
}
\examples{

load(system.file('extdata', 'data.RData', package='DEsubs'))

outfile <- tempfile(fileext='.pdf')  

res <- subpathwayToGraph( 
                DEsubs.out=DEsubs.out, 
                submethod='community.walktrap', 
                subname=paste0('sub', 6),
                size=c(10,10),
                export='pdf',
                outfile=outfile )

}
