

\name{subpathwayVisualization}
\alias{subpathwayVisualization}
\title{Subpathway level visualization}
\usage{
subpathwayVisualization( DEsubs.out, references, submethod, subname, colors, 
                        scale, shuffleColors, outfiles, export, verbose )
}
\arguments{
\item{DEsubs.out}{Return value from \code{\link{DEsubs}}.}
\item{references}{Topological references include degree, betweeness centrality,
    closeness centrality, hub_score, eccentricity and page_rank. 
    The corresponding options are 'degree', 'betweenness',  'closeness', 
    'hub_score', 'eccentricity' and 'page_rank'. Functional references include
    (i) KEGG's pathway terms, (ii) Gene Ontologies of Molecular function, 
    biological processes and cellular components) (iii) Disease terms 
    from OMIM and GAD databases), (iv) Drug substances from DrugBank) and 
    the influence of (v)  microRNA targets from miRecords and (vi) 
    Transcription Factor targets from Transfac and Jaspar. The corresponding 
    options are 'KEGG', 'GO_bp',  'GO_cc', 'GO_mf', 'Disease_OMIM' ,
    'Disease_GAD', 'Drug_DrugBank', 'miRNA' and 'TF'.Using option 'all' 
    results in the selection of all afforementioned options, along with any 
    other custom gene sets within the 'DEsubs/Data' user specified-directory.}
\item{submethod}{Subpathway extraction type selection 
    (see all 124 options along with their R commands in supplementary 
    document)}
\item{subname}{Subpathway name as contained in \code{\link{DEsubs}} return 
value, i.e. 'sub1'}
\item{colors}{ A custom color mode which overrrides the default settings.}
\item{scale}{ A value in (0,1] used to scale the visualization. Useful in the 
case of long labels which are trimmed by default.}
\item{shuffleColors}{ TRUE to shuffle user defined or default colors. 
Defaults to FALSE. }
\item{outfiles}{ Output filenames of the visualizations. If the argument is 
not specified, default filenames are used ('DEsubs/Output').}
\item{export}{ Export type of visualizations ('plot', 'pdf') }
\item{verbose}{ TRUE to display informative messages, FALSE to hide. }
}
\value{
A list of matrices, each containing the  P-Value of enrichment between the 
terms for a specific reference (rows) and each of the subpathway genes 
(columns). If there is no enrichment, the value is NA. 
}
\description{
Circular diagrams containing subpathways enrichment in potential key 
regulators (miRNAs, TFs) and biological, biomedical and pharmacological 
issues.
}
\details{The associations of subpathways with various biological and 
    pharmacologicalfeatures are estimated through a hypergeometric test. 
    The enriched associations of a subpathway to each feature are 
    illustrated through circular diagrams.
}
\examples{

load(system.file('extdata', 'data.RData', package='DEsubs'))

outfile <- tempfile(fileext='.pdf') 

res <- subpathwayVisualization( 
                DEsubs.out=DEsubs.out,  
                references=c('TF'), 
                submethod='community.walktrap',
                subname='sub1', 
                scale=c(1.0), 
                export='pdf',
                outfile=outfile )

}

