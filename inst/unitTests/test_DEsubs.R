test_DEsubs <- function() 
{
   message('Testing DEsubs...', appendLF=FALSE)

	load(system.file('extdata', 'data.RData', package='DEsubs'))

	DEsubs.run <- DEsubs(   org='hsa', 
	                        mRNAexpr=mRNAexpr, 
	                        mRNAnomenclature='entrezgene', 
	                        pathways='All', 
	                        DEtool=NULL, 
	                        DEpar=0.05,
	                        CORtool='pearson', 
	                        CORpar=0.6, 
	                        subpathwayType=NULL,
	                        rankedList=rankedList)

   all.equal(DEsubs.run, DEsubs.out)

   message('done.')
}
