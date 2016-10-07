

\name{organismVisualization}
\alias{organismVisualization}
\title{Organism level visualization}
\usage{
organismVisualization( DEsubs.out, references, topSubs, topTerms, colors,
                        export, width, height, outfiles, verbose )
}
\arguments{
\item{DEsubs.out}{Return value from \code{\link{DEsubs}} }
\item{references}{Functional associations with 
    (i) KEGG's pathway terms, (ii) Gene Ontologies of Molecular function, 
    biological processes and cellular components) (iii) Disease terms 
    from OMIM and GAD databases), (iv) Drug substances from DrugBank) and 
    the influence of (v)  microRNA targets from miRecords and (vi) 
    Transcription Factor targets from Transfac and Jaspar. The corresponding 
    options are 'KEGG', 'GO_bp',  'GO_cc', 'GO_mf', 'Disease_OMIM' ,
    'Disease_GAD', 'Drug_DrugBank', 'miRNA' and 'TF'. If no argument is 
    supplied, no option is selected. Using option 'all' results in the 
    selection of all afforementioned options, along with any other custom
    gene sets within the 'DEsubs/Data' user specified-directory.}
\item{topSubs}{ Default value is 10 }
\item{topTerms}{ Default value is 20 }
\item{colors}{ A custom color mode which overrides the default settings. }
\item{export}{ Export type of visualizations ('plot', 'pdf') }
\item{width}{ The width of the printable area (pdf) }
\item{height}{ The height of the printable area (pdf) }
\item{outfiles}{ Output filenames of the visualizations. If the argument is 
not specified, default filenames are used ('DEsubs/Output').}
\item{verbose}{ TRUE to display informative messages, FALSE to hide. }
}
\value{
No value is returned.
}
\description{
Organism level measures
}


\examples{

load(system.file('extdata', 'data.RData', package='DEsubs'))

res <- organismVisualization( 
                DEsubs.out=DEsubs.out, references='KEGG', 
                topSubs=10, topTerms=20,
                export='pdf', width=7, height=6,
                verbose=FALSE)

}
