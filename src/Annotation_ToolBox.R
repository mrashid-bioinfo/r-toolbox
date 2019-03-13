## ================================================================== ##
##		Extract Gene Chr, Start and End from Ensembl BioMart
##
## ================================================================== ##


gene_annotation_biomart = function( biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org", dataset = "hsapiens_gene_ensembl", genes = c("BRCA1", "TP53") )
{
	library(biomaRt)
	ens37mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org", path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl") 
	filters <- listFilters(ens37mart)
	attributes <- listAttributes(ens37mart)
	rec_genes.coords37 <- getBM(attributes=c( "hgnc_symbol",  "start_position", "end_position", "chromosome_name"), filters  = 'hgnc_symbol',  values = genes, mart = ens37mart) 
	rec_genes.coords37 = rec_genes.coords37[ -grep("HG",rec_genes.coords37$chromosome_name), ]
	#gene_locs = paste(rec_genes.coords37[,4], ":", rec_genes.coords37[,2], "-", rec_genes.coords37[,3], sep = "" )

	return(gene_locs)
}

## ------------------------------- ##
##
## 	Function :: gene_exon_biomart
## 	Returns gene exons boundaries
##	
##
##
gene_exon_biomart = function( biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org", dataset = "hsapiens_gene_ensembl", genes = c("BRCA1", "TP53") )
{
	library(biomaRt)
	ens37mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org", path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl") 
	filters <- listFilters(ens37mart)
	attributes <- listAttributes(ens37mart)
	rec_genes.coords37 <- getBM(attributes=c( "external_gene_name", "ensembl_exon_id", "rank", "chromosome_name",  "exon_chrom_start", "exon_chrom_end"), filters  = 'hgnc_symbol',  values = genes, mart = ens37mart) 
	rec_genes.coords37 = rec_genes.coords37[ grep("HG",rec_genes.coords37$chromosome_name, invert = T ), ]
	rec_genes.coords37_uniq = rec_genes.coords37[ !duplicated( rec_genes.coords37$rank ), ]
	#gene_locs = paste(rec_genes.coords37[,4], ":", rec_genes.coords37[,2], "-", rec_genes.coords37[,3], sep = "" )

	return(rec_genes.coords37_uniq)
}

cBio_conseq_mapping = function( conseq = NULL )
{
	if( conseq == "missense" ) return( "Missense_Mutation" )
	if( conseq == "nonsense" ) return( "Nonsense_Mutation" )
	if( conseq == "ess_splice" ) return( "Splice_Site" )
}

cBio_splice_mapping = function( chr = NULL, loc = NULL )
{
	require(biomaRt)
	ens37mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org", path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl") 
	filters <- listFilters(ens37mart)
	attributes <- listAttributes(ens37mart)
	
	gene_name = apply( query_df, 1, function(X){
		temp_list = list(chromosome_name = X[1], start = X[2], end = X[3] )
		gene = getBM(attributes=c( "hgnc_symbol"), filters  = c('chromosome_name', 'start', 'end'), values = temp_list, mart = ens37mart, uniqueRows = TRUE )
		gene_collapsed = paste( gene[,1], collapse = "," );
	} )
}

mutation_mapper_format = function( df_file = NULL, gene_name_colname = "Gene_Name", gene_name = "ALPK1", sample_tag = "PD", conseq_colname = "Deletorious_Consequence", conseq_summary_colname = "Consequence_Summary" , variant_type_colname = "Variant_Type")
{
	require(data.table);
	
	df = read.delim( df_file, sep = "\t", header = T, stringsAsFactors = F, check.names = F );
	temp_df = df[ df[,gene_name_colname] == gene_name, ];

	sampleIdx = grep( paste("^", sample_tag, sep = ""), colnames(temp_df) )
	sampleIDs = colnames(temp_df)[ sampleIdx ]
	conseq_colIdx = grep( conseq_colname, colnames(temp_df) )

	cBio_df = c()
	for( i in 1:dim(temp_df)[1])
	{
		x = as.numeric(temp_df[i, sampleIDs])
		mutated_samples = paste( sampleIDs [ x != 0 ], collapse = "," )
		protein_change = unlist ( strsplit( temp_df[i, conseq_colIdx], split = "|", fixed  = T ) )[5]
		protein_change = gsub("p.", "", protein_change )
		conseq = cBio_conseq_mapping( temp_df[i, conseq_summary_colname] )
		chr = temp_df[i, 1 ] 
		start = temp_df[i, 2 ] 
		end = temp_df[i, 3 ] 
		ref = temp_df[i, 4 ] 
		alt = temp_df[i, 5 ] 

		## Take the closest exon protein change ##
		if( temp_df[i, conseq_summary_colname] ==  "ess_splice" )
		{
			protein_change = "Splice_Site"
			conseq = "Splice_Site"
		}

		if( temp_df[i, conseq_summary_colname] ==  "Frameshift" )
		{
			if( temp_df[i, variant_type_colname] == "Del"  ){ conseq = "Frame_Shift_Del" }
			if( temp_df[i, variant_type_colname] == "Ins"  ){ conseq = "Frame_Shift_Ins" }
			
			protein_change = gsub( "\\*[0-9]+", "", protein_change);
		}

		if( conseq %in% c( "Missense_Mutation", "Nonsense_Mutation", "Splice_Site" ) )
		{
			start = end
		}

		cBio_string = c( gene_name, mutated_samples, protein_change, conseq, chr, start, end, ref, alt );
		cBio_df = rbind( cBio_df, cBio_string )	
	}

	cBio_df = as.data.frame(cBio_df);
	colnames(cBio_df) = c( "Hugo_Symbol","Sample_ID","Protein_Change","Mutation_Type","Chromosome","Start_Position","End_Position","Reference_Allele","Variant_Allele" )
	cBio_df_E = cBio_df %>% mutate( Sample_ID = strsplit( as.character(Sample_ID), split = ",") ) %>% unnest(Sample_ID) 

	cBio_df_F = data.frame( Hugo_Symbol = cBio_df_E[,"Hugo_Symbol"], Sample_ID = cBio_df_E[,"Sample_ID"], cBio_df_E[,c("Protein_Change","Mutation_Type","Chromosome","Start_Position","End_Position","Reference_Allele","Variant_Allele")] )


	return( cBio_df_F );

} ## End of function ##


