

\name{DEsubs}
\alias{DEsubs}
\title{Default run of DEsubs}
\usage{
DEsubs(org, mRNAexpr, mRNAnomenclature, pathways, DEtool, DEpar,
    CORtool, CORpar, subpathwayType, rankedList, verbose)
}
\arguments{
\item{org}{Organism identifier ('hsa', 'mmu', 'rno', 'sce', 'ath', 'dme', 
'dre')}
\item{mRNAexpr}{RNA-seq expression data in the form of either a matrix
or a filename of a text file stored in the 'User' directory.}
\item{mRNAnomenclature}{mRNAnomenclature ('entrezgene', 'ensembl_gene_id',
    'ensembl_transcript_id', 'ensembl_peptide_id', 'hgnc_id', 'hgnc_symbol',
    'hgnc_transcript_name', 'refseq_mrna', 'refseq_peptide')}
\item{pathways}{Pathway type ('All', 'Non-Metabolic', 'Metabolic')}
\item{DEtool}{DEG analysis tool selection for NodeRule ('edgeR', 'DESeq', 
    'EBSeq', 'samr', 'NBPSeq', 'voom+limma', 'vst+limma', 'TSPM')}
\item{DEpar}{DE analysis tools Q-value threshold of NodeRule 
    (default: DEGpar = 0.05)}
\item{CORtool}{ Correlation measure selection for EdgeRule 
    ('pearson', 'kendall', 'spearman') }
\item{CORpar}{ Correlation measure threshold of EdgeRule 
            (default: CORpar = 0.6) }
\item{subpathwayType}{ Subpathway extraction type selection 
    (get all availiable options from \link{subpathwayTypes}) }
\item{rankedList}{ A named vector of genes and their corresponding 
significance of differential expression in the form of a Q-value. If the 
argument is not null, no DEtool is used for differential expression analysis.}
\item{verbose}{ TRUE to display informative messages, FALSE to hide. }
}
\value{
A list used as input in  
    \code{\link{geneVisualization}}, \code{\link{subpathwayVisualization}}, 
    \code{\link{subpathwayToGraph}}, \code{\link{organismVisualization}}
}
\description{
Default run of DEsubs 
}
\details{
\itemize{
\item{ Class vector needed values 1 and 2 for each class 
        (eg. control and disease samples). }
\item{ DEpar should be less than 0.05 in order to return statistically 
        significant DEGs.}
\item{ Higher CORpar values result in stricter correlation criteria, 
        i.e. less acceptable interactions. }
}
}
\examples{

load(system.file('extdata', 'data.RData', package='DEsubs'))

DEsubs.run <- DEsubs(   org='hsa', 
                        mRNAexpr=mRNAexpr, 
                        mRNAnomenclature='entrezgene', 
                        pathways='All', 
                        DEtool=NULL, DEpar=0.05,
                        CORtool='pearson', CORpar=0.6, 
                        subpathwayType=NULL,
                        rankedList=rankedList)

}
