	genes = c("VHL","IDH1","RB1","PIK3CA");
	library(biomaRt)
	ens37mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org", path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl") 
	filters <- listFilters(ens37mart)
	attributes <- listAttributes(ens37mart)
	
	df = c()
	for( i in 1:length(genes) )
	{
	
		## Exon data
		genes.coords37 <- getBM( attributes = c( "chromosome_name", "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end", "rank" ), filters  = 'hgnc_symbol', values = genes[i], mart = ens37mart );
		genes.coords37_uniq = genes.coords37 [ !duplicated(genes.coords37$rank), ]

		df = rbind( df, genes.coords37_uniq )
	}