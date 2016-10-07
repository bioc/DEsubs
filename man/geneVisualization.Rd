

\name{geneVisualization}
\alias{geneVisualization}
\title{Gene level visualization}
\usage{
geneVisualization(DEsubs.out, measures.topological, measures.functional, 
    measures.barplot, topGenes, colors.topological, colors.functional, 
    colors.barplot, size.topological, size.functional, size.barplot, 
    outfile.topological, outfile.functional, outfile.barplot,
    export, verbose)
}
\arguments{
\item{DEsubs.out}{Return value from \code{\link{DEsubs}} }
\item{measures.topological}{Functional visualization type(s). }
\item{measures.functional}{Topological visualization type(s). }
\item{measures.barplot}{Gene level Visualization type}
\item{topGenes}{Number of genes with greater Q-values. Default value is 10.  }
\item{colors.topological}{ A custom color mode which overrrides the 
default settings.}
\item{colors.functional}{ A custom color mode which overrrides the default 
settings.}
\item{colors.barplot}{ A custom color mode which overrrides the default 
settings.}
\item{size.topological}{ A vector storing width and height of the topological
measures visualization}
\item{size.functional}{ A vector storing width and height of the functional
measures visualization }
\item{size.barplot}{ A vector storing width and height of the barplot
visualization }
\item{outfile.topological}{ Output file name of the topological measures
visualization. }
\item{outfile.functional}{ Output file name of the functional measures
visualization. }
\item{outfile.barplot}{ Output file name of the barplot visualization. }
\item{export}{ Export type of visualizations ('plot', 'pdf') }
\item{verbose}{ TRUE to display informative messages, FALSE to hide. }
}
\value{
Individual measure results in matrix form.
}
\description{
Visualizes topologically and functionally significant genes using graph 
theory measures as well as their correlation to pathway, disease, drug, 
ontology, microRNA and Transcription Factor terms based on external 
references.
}
\details{
\itemize{
\item{ Topological visualization type contains six topological graph theory 
    measures, namely degree, betweeness centrality, closeness centrality, 
    hub_score, eccentricity and page_rank. The availiable options are
    'degree', 'betweenness',  'closeness', 'hub_score', 'eccentricity' and
    'page_rank'. If no argument is supplied, all availiable options are 
    selected by default. If the argument is NULL, no visualization is 
    exported.}
\item{ Functional visualization type contains associations with 
    (i) KEGG's pathway terms, (ii) Gene Ontologies of Molecular function, 
    biological processes and cellular components) (iii) Disease terms 
    from OMIM and GAD databases), (iv) Drug substances from DrugBank) and 
    the influence of (v)  microRNA targets from miRecords and (vi) 
    Transcription Factor targets from Transfac and Jaspar. The availiable 
    options are 'KEGG', 'GO_bp',  'GO_cc', 'GO_mf', 'Disease_OMIM' ,
    'Disease_GAD', 'Drug_DrugBank', 'miRNA' and 'TF'. If no argument is 
    supplied, all availiable options are selected by default. If the argument
    is NULL, no visualization is exported. Using option 'all' results in the 
    selection of all afforementioned options, along with any other custom
    gene sets within the 'DEsubs/Data' user specified-directory.}
}
}
\examples{

load(system.file('extdata', 'data.RData', package='DEsubs'))

res <- geneVisualization(  
            DEsubs.out=DEsubs.out, top=10,
            measures.topological=c( 'degree', 'betweenness', 'closeness',
                                    'hub_score', 'eccentricity', 'page_rank'),
            measures.functional=c(  'KEGG', 
                                    'Disease_OMIM', 'Disease_GAD', 
                                    'Drug_DrugBank','miRNA', 'TF'),
            size.topological=c(5,4),
            size.functional=c(7,4),
            size.barplot=c(5,6),
            export='pdf', verbose=FALSE)

}
