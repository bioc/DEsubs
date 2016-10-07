

\name{subpathwayTypes}
\alias{subpathwayTypes}
\title{All subpathway types}
\usage{
subpathwayTypes(grouping)
}
\arguments{
\item{grouping}{By supplying one of the availiable groupings, specific 
subsets of availiable subpathway types can be extracted.}
}
\value{
A vector of all 124 basic subpathway types. See Details section for handly 
groupings.
}
\description{
All subpathway types.
}
\details{
Apart from the 124 distinct subpathway types, several groupings are availiable:
\describe{
\item{'all'}{ All subpathway options. }
\item{'all.bwd'}{ Forward propagation. }
\item{'all.fwd'}{Backward propagation. }
\item{'all.stream'}{All possible genes targeting or are targeted 
via a path starting from a gene of interest.}
\item{'all.neighbourhood'}{All adjacent genes of a gene of interest with 
incoming or outgoing links.}
\item{'all.cascade'}{A finite sequence of interactions connecting a sequence 
of genes starting or ending from a gene of interest.}
\item{'all.community'}{Group of genes sharing common properties.}
\item{'all.component'}{Subgraphs on which any two vertices are (strongly)
connected to each other by paths.}
\item{ }{}
\item{'all.functional'}{Forward and backward streams starting from 
genes/nodes with crucial functional role within the network. Individual 
options include GO_bp, GO_cc, GO_mf, KEGG, Disease_OMIM, Disease_GAD, 
Drug_DrugBank, miRNA, TF, DEG, which are defined below.}
\item{'all.GO_bp'}{Genes acting as a bridge among Gene Ontology (GO) 
Biological Process terms. }
\item{'all.GO_cc'}{Genes acting as a bridge among Gene Ontology (GO) 
Cellular Component terms. }
\item{'all.GO_mf'}{Genes acting as a bridge among Gene Ontology (GO) 
Molecular Function terms. }
\item{'all.KEGG'}{Genes acting as a bridge among KEGG pathway maps. }
\item{'all.Disease_OMIM'}{Genes acting as a bridge among OMIM Disease 
targets.}
\item{'all.Disease_GAD'}{Genes acting as a bridge among GAD Disease targets.}
\item{'all.Drug_DrugBank'}{Genes acting as a bridge among DrugBank Drug 
targets.}
\item{'all.miRNA'}{Genes acting as a bridge among microRNA-gene targets. }
\item{'all.TF'}{Genes acting as a bridge among TF-gene targets. }
\item{'all.DEG'}{Genes with highly Differentially expressed by each 
experimental data. }
\item{}{}
\item{'all.topological'}{Forward and backward streams starting from 
genes/nodes with crucial topological role within the network. Individual 
options include degree, betweenness, closeness, hub_score, eccentricity, 
page_rank, start_nodes which are defined below.}
\item{'all.degree'}{Number of gene's adjacent interactions.}
\item{'all.betweenness'}{Number of shortest paths from all vertices to  
others passing through a node.}
\item{'all.closeness'}{Inverse of farness, which is the sum of distances 
to all other nodes.}
\item{'all.hub_score'}{Kleinberg's hub centrality score.}
\item{'all.eccentricity'}{Shortest path distance from the farthest other 
node in the graph.}
\item{'all.page_rank'}{Google Page Rank.}
\item{'all.start_nodes'}{Gene nodes without incoming links.}
\item{}{}
}

An exhaustive list of all 124 subpathway types follows:

\tabular{ll}{
    \tab  \cr
    STREAM-TOPOLOGICAL  \tab  \cr
    'bwd.stream.topological.degree' \tab  
    'fwd.stream.topological.degree' \cr 
    'bwd.stream.topological.betweenness' \tab 
    'fwd.stream.topological.betweenness'  \cr
    'bwd.stream.topological.closeness' \tab 
    'fwd.stream.topological.closeness' \cr 
    'bwd.stream.topological.hub_score' \tab 
    'fwd.stream.topological.hub_score' \cr
    'bwd.stream.topological.eccentricity' \tab 
    'fwd.stream.topological.eccentricity' \cr 
    'bwd.stream.topological.page_rank' \tab  
    'fwd.stream.topological.page_rank' \cr
    'bwd.stream.topological.start_nodes' \tab 
    'fwd.stream.topological.start_nodes' \cr 
    \tab  \cr
    STREAM-FUNCTIONAL  \tab  \cr
    'bwd.stream.functional.GO_bp' \tab 
    'fwd.stream.functional.GO_bp' \cr
    'bwd.stream.functional.GO_cc' \tab 
    'fwd.stream.functional.GO_cc' \cr 
    'bwd.stream.functional.GO_mf' \tab 
    'fwd.stream.functional.GO_mf' \cr
    'bwd.stream.functional.KEGG' \tab 
    'fwd.stream.functional.KEGG' \cr 
    'bwd.stream.functional.Disease_OMIM' \tab 
    'fwd.stream.functional.Disease_OMIM' \cr
    'bwd.stream.functional.Disease_GAD' \tab 
    'fwd.stream.functional.Disease_GAD' \cr
    'bwd.stream.functional.Drug_DrugBank' \tab 
    'fwd.stream.functional.Drug_DrugBank' \cr
    'bwd.stream.functional.miRNA' \tab 
    'fwd.stream.functional.miRNA' \cr
    'bwd.stream.functional.TF' \tab 
    'fwd.stream.functional.TF' \cr
    'bwd.stream.functional.DEG' \tab 
    'fwd.stream.functional.DEG' \cr
    \tab  \cr
    NEIGHBOURHOOD-TOPOLOGICAL  \tab  \cr
    'bwd.neighbourhood.topological.degree' \tab 
    'fwd.neighbourhood.topological.degree' \cr
    'bwd.neighbourhood.topological.betweenness' \tab 
    'fwd.neighbourhood.topological.betweenness' \cr
    'bwd.neighbourhood.topological.closeness' \tab 
    'fwd.neighbourhood.topological.closeness' \cr
    'bwd.neighbourhood.topological.hub_score' \tab 
    'fwd.neighbourhood.topological.hub_score' \cr
    'bwd.neighbourhood.topological.eccentricity' \tab 
    'fwd.neighbourhood.topological.eccentricity' \cr
    'bwd.neighbourhood.topological.page_rank' \tab 
    'fwd.neighbourhood.topological.page_rank' \cr
    'bwd.neighbourhood.topological.start_nodes' \tab 
    'fwd.neighbourhood.topological.start_nodes' \cr
    \tab  \cr
    NEIGHBOURHOOD-FUNCTIONAL   \tab  \cr
    'bwd.neighbourhood.functional.GO_bp' \tab 
    'fwd.neighbourhood.functional.GO_bp' \cr
    'bwd.neighbourhood.functional.GO_cc' \tab 
    'fwd.neighbourhood.functional.GO_cc' \cr
    'bwd.neighbourhood.functional.GO_mf' \tab 
    'fwd.neighbourhood.functional.GO_mf' \cr
    'bwd.neighbourhood.functional.KEGG' \tab 
    'fwd.neighbourhood.functional.KEGG' \cr
    'bwd.neighbourhood.functional.Disease_OMIM' \tab 
    'fwd.neighbourhood.functional.Disease_OMIM' \cr
    'bwd.neighbourhood.functional.Disease_GAD' \tab 
    'fwd.neighbourhood.functional.Disease_GAD' \cr
    'bwd.neighbourhood.functional.Drug_DrugBank' \tab 
    'fwd.neighbourhood.functional.Drug_DrugBank' \cr
    'bwd.neighbourhood.functional.miRNA' \tab 
    'fwd.neighbourhood.functional.miRNA' \cr
    'bwd.neighbourhood.functional.TF' \tab 
    'fwd.neighbourhood.functional.DEG' \cr
    \tab  \cr
    CASCADE-TOPOLOGICAL  \tab  \cr
    'bwd.cascade.topological.degree' \tab 
    'fwd.cascade.topological.degree' \cr
    'bwd.cascade.topological.betweenness' \tab 
    'fwd.cascade.topological.betweenness' \cr
    'bwd.cascade.topological.closeness' \tab 
    'fwd.cascade.topological.closeness' \cr
    'bwd.cascade.topological.hub_score' \tab 
    'fwd.cascade.topological.hub_score' \cr
    'bwd.cascade.topological.eccentricity' \tab 
    'fwd.cascade.topological.eccentricity' \cr
    'bwd.cascade.topological.page_rank' \tab 
    'fwd.cascade.topological.page_rank' \cr
    'bwd.cascade.topological.start_nodes' \tab 
    'fwd.cascade.topological.start_nodes' \cr
    \tab  \cr
    \tab  \cr
    CASCADE-FUNCTIONAL  \tab  \cr
    'bwd.cascade.functional.GO_bp' \tab 
    'fwd.cascade.functional.GO_bp' \cr
    'bwd.cascade.functional.GO_cc' \tab 
    'fwd.cascade.functional.GO_cc' \cr
    'bwd.cascade.functional.GO_mf' \tab 
    'fwd.cascade.functional.GO_mf' \cr
    'bwd.cascade.functional.KEGG' \tab 
    'fwd.cascade.functional.KEGG' \cr
    'bwd.cascade.functional.Disease_OMIM' \tab 
    'fwd.cascade.functional.Disease_OMIM' \cr
    'bwd.cascade.functional.Disease_GAD' \tab 
    'fwd.cascade.functional.Disease_GAD' \cr
    'bwd.cascade.functional.Drug_DrugBank' \tab 
    'fwd.cascade.functional.Drug_DrugBank' \cr
    'bwd.cascade.functional.miRNA' \tab 'fwd.cascade.functional.miRNA' \cr
    'bwd.cascade.functional.TF' \tab 'fwd.cascade.functional.TF' \cr
    'bwd.cascade.functional.DEG' \tab 'fwd.cascade.functional.DEG' \cr
    \tab  \cr
    COMMUNITY  \tab  \cr
    'community.walktrap' \tab 'community.edge_betweenness' \cr
    'community.fast_greedy' \tab 'community.leading_eigen' \cr
    'community.infomap' \tab 'community.louvain' \cr
    \tab  \cr
    COMPONENT-CLIQUES  \tab  \cr
    'component.max_cliques' \tab 'component.decompose' \cr
    'component.3-cliques' \tab ... \cr
    ... \tab 'component.9-cliques' \cr
    \tab  \cr
    COMPONENT-CORENESS  \tab  \cr
    'component.3-coreness' \tab ... \cr
    ... \tab 'component.9-coreness' \cr
}

}
\examples{

basic.types <- subpathwayTypes()
stream.types <- subpathwayTypes(grouping='all.stream')

}








