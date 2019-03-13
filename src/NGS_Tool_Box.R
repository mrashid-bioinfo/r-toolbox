 ###################################
 #                                 #
 ##  NGS Data Process R Library  ##
 #                                 #
 ###################################

###################################
##
##	Convert a genotype dataframe to binary data frame
##
###################################

genotype_binary_conversion = function( genotype, gt_ref = "0/0", gt_nocall = "./." )
{
	return( ifelse( genotype == gt_ref | genotype == gt_nocall, 0, 1 ) )
}



###################################
##
##	Find the most severe consequence and gene from VEP csq tag 
##
##
###################################

merge_snv_and_indel = function( genes_snv_df = genes_snv_df, genes_indel_df = genes_indel_df, sample_tag = "^PD" )
{
	## Merge SNV and Indel data
	snv_samples = colnames(genes_snv_df)[ grep(sample_tag, colnames(genes_snv_df)) ]
	indel_samples = colnames(genes_indel_df)[ grep(sample_tag, colnames(genes_indel_df)) ]
	common_samples = intersect( snv_samples, indel_samples )
	uniq_genes = unique( c( genes_indel_df$Gene_Name, genes_snv_df$Gene_Name ) )
	snv_indel_gene = matrix(0, length(uniq_genes), length(common_samples) + 1 )
	colnames( snv_indel_gene) = c( "Count",  common_samples );
	rownames(snv_indel_gene) = uniq_genes;
	
	for( i in 1:length(uniq_genes) )
	{
		print(uniq_genes[i])
		gene = uniq_genes[i];
		for( j in 1:length(common_samples) )
		{
			print(common_samples[j])
			sample = common_samples[j];
			snv_val = genes_snv_df[ gene, sample ];
			indel_val = genes_indel_df[ gene, sample ];
			if( as.character(indel_val) == "0" || is.na(indel_val) )
			{
				new_val = as.character(snv_val)	
			} 
			else if( as.character(snv_val) == "0" && as.character(indel_val) == "0" )
			{
				new_val = "0"
			}
			else if( as.character(snv_val) == "0" || is.na(snv_val) )
			{
				new_val = as.character(indel_val)	
			} 
			else{
				new_val = c( as.character(snv_val), as.character(indel_val) )
			}
			snv_indel_gene[ gene, sample] =  paste( new_val , collapse = "," )
		}
		snv_indel_gene[gene, "Count"] = length( snv_indel_gene[ gene, ] [ snv_indel_gene[ gene, ] != "0" ] )
	}

	write.table(snv_indel_gene, "snv_indel_gene.txt", sep = "\t", quote = F)

	## Fix : Replace NA to 0

	snv_indel_gene_mod1 = t( apply( snv_indel_gene, 1, function(X)
	{
		ifelse( X=="NA", 0, X )
	} ) )
	write.table(snv_indel_gene_mod1, "snv_indel_gene_mod1.txt", sep = "\t", quote = F)

	snv_indel_gene[,"Count"] = apply(snv_indel_gene[,common_samples], 1, function( X ){ length( X[ X != "NA" & X != "0"] ) } )
	snv_indel_gene_o = snv_indel_gene[ order( as.numeric( snv_indel_gene[, "Count"] ), decreasing = T), ]

	## With more than 2 sample recurrence 
	snv_indel_gene_rec = snv_indel_gene_o[ as.numeric(snv_indel_gene_o[,"Count"]) >2,  ]

	return( snv_indel_gene_o )

}

###################################
##
##	Find the most severe consequence and gene from VEP csq tag 
##
##	input : VEP CSQ tag
##	output : [ Consequence and gene]
##

severe_conseq = function( csq = NULL )
{
	csq_split = unlist ( strsplit( csq, split = ",", fixed = F ) );
	if( length( csq_split[ grep( "HIGH", csq_split) ] ) > 0 )
	{
		severe = csq_split[ grep( "HIGH", csq_split)[1] ]	
	}
	else if( length( csq_split[ grep( "MODERATE", csq_split) ] ) > 0 )
	{
		severe = csq_split[ grep( "MODERATE", csq_split)[1] ]	
	}
	else if( length( csq_split[ grep( "LOW", csq_split) ] ) > 0 )
	{
		severe = csq_split[ grep( "LOW", csq_split)[1] ]	
	}
	else if( length( csq_split[ grep( "MODIFIER", csq_split) ] ) > 0 )
	{
		severe = csq_split[ grep( "MODIFIER", csq_split)[1] ]	
	}

	print(severe)
	temp = unlist( strsplit( severe, split = "|", fixed = T ) );

	return( c( temp[2], temp[3], temp[4] ) )
}



## --------------------------------------------------------------- ##
## Function to convert VCF to Dataframe and select non syn mutations
##	
##	Takes a bedR split data frame. Find the most deletorious conseq
##
vcf_to_conseq_split_df = function( df )
{
	conseq = sapply(df$vcf$CSQ, function(X)
	{ 
		a = unlist( strsplit( X, split = "," , fixed = T ) );  
		temp_conseq = c()
		for( i in 1:length( a ) )
		{
			b = as.character ( unlist ( strsplit( a[i], split = "|", fixed = T ) ) )[2]
			temp_conseq = append( temp_conseq, b )
		}
		return(paste( temp_conseq, collapse = "," ))
	} )

	vcf_df = data.frame( df$vcf, Consequence = conseq,  stringsAsFactors = F )
	vcf_df_E = vcf_df %>% mutate( Consequence= strsplit( Consequence, split = "," ) ) %>% unnest( Consequence )
	vcf_df_E$id = return_id( vcf_df_E, c(1,2,4,5,49) )
	vcf_df_E_uniq = vcf_df_E[ !duplicated( vcf_df_E$id ), ]
	vcf_df_mod1 = separate( vcf_df_E_uniq, col = "id", into = c("Chr1","Pos1","Ref1","Alt1","Consequence"), sep = "-" )
	vcf_df_mod2 = vcf_df_mod1[ , c(1:48,53) ]
	vcf_df_mod2$Consequence = variant_remap(vcf_df_mod2$Consequence)
	vcf_df_mod2_nonsyn = vcf_df_mod2[ vcf_df_mod2$Consequence %in% conseq_to_keep , ]
	#vcf_df_mod2_nonsyn_id = return_id( vcf_df_mod2_nonsyn, c(1,2,4,5), collapse = "_" )
	
	return( vcf_df_mod2_nonsyn )

}


## ==================================================================== ##
##
## 	Function :: Tumour normal read count summary   						##
##	
##	Given two CGP format read count metrices it returns a list
##	containing two metrices summing all read count for all possible 
##	bases. 
##	
##	Input :: 
##		1. Two metrices
##		
##	output :: One list containing two metrices
##

get_tumour_normal_read_count = function( tumour_depth_df = NULL, normal_depth_df = NULL )
{
	tumour_sampleIds = colnames(tumour_depth_df)[grep("^PD", colnames( tumour_depth_df)) ];
	normal_sampleIds = colnames(normal_depth_df)[grep("^PD", colnames( normal_depth_df)) ];

	tumour_read_mat = matrix( NA, dim(tumour_depth_df)[1], length(tumour_sampleIds) )
	normal_read_mat = matrix( NA, dim(tumour_depth_df)[1], length(normal_sampleIds) )

	for( i in 1:length(tumour_sampleIds) )
	{
		temp_tumour = tumour_depth_df[,tumour_sampleIds[i]]
		tumour_max_count = sapply( temp_tumour, function(X) { sum( as.numeric( unlist( strsplit( X, split= "|", fixed = T ) ) ) )  } ); 
		tumour_read_mat[ , i ] = tumour_max_count;

		temp_normal = normal_depth_df[,normal_sampleIds[i]]
		normal_max_count = sapply( temp_normal, function(X) { sum( as.numeric( unlist( strsplit( X, split= "|", fixed = T ) ) ) )  } ); 
		normal_read_mat[ , i ] = normal_max_count;
	}

	colnames(tumour_read_mat) = tumour_sampleIds
	colnames(normal_read_mat) = normal_sampleIds

	mat_list = list( tumour = tumour_read_mat, normal = normal_read_mat );
	return( mat_list );
}



## -------------------------------------- ##
##
## Function :: get_mutated_samples
##
##	Parameters :: 
##
##
## For every SNV get the mutated samples from a CGP formated TAB delimited file 
##
get_mutated_samples = function( df = NULL, val = 0, tag = "PD" )
{
	sampleIDs = colnames(df)[ grep( tag, colnames(df) ) ]
	mutated_samples = apply( df[,sampleIDs], 1, function(X) { return(paste(sampleIDs [ X != val], collapse = "," )) } )
	return(mutated_samples)
}



germline_check = function( normal_ref_read, normal_alt_read, tumour_ref_read, tumour_alt_read )
{

	norm_total = as.numeric(normal_alt_read) + as.numeric(normal_ref_read )
	if(  norm_total > 0 )
	{
		normal_alt_frac = as.numeric(normal_alt_read) / norm_total;	
	}
	else{
		normal_alt_frac = 0
	}
	

	tum_total = as.numeric(tumour_alt_read) + as.numeric(tumour_ref_read )
	if( tum_total > 0 ) 
	{
		tumour_alt_frac = as.numeric(tumour_alt_read) / tum_total;	
	}
	else{
		tumour_alt_frac = 0
	}

	if( as.numeric( tumour_alt_read ) == 0 )
	{
		verdict = "No Evidence"	
	}
	else if( as.numeric( tumour_alt_read ) <= 10 )
	{
		verdict = "Low Evidence"
	}
	else if( as.numeric( normal_alt_frac ) >= as.numeric( tumour_alt_frac ) )
	{
		verdict = "Germline"
	}
	else if( as.numeric( normal_alt_frac ) >= 0.01 )
	{
		verdict = "Germline"
	}
	else if( as.numeric( tumour_alt_frac ) >= 0.10 )
	{
		verdict = "Somatic"
	}
	else{
		verdict = "Not Sure"	
	}

	return (verdict)
}

## Maps variant consequences to known terms ##
## Input :: Character Vector [ e.g. c("missense_variant","stop_gained" ... ) ]
##
variant_remap = function( variant_conseq )
{
	variant_conseq_remapped = c()
	for( i in 1:length( variant_conseq ) )
	{	
		if(  length( grep( "stop_gained", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Nonsense"	
		}
		else if(  length( grep( "stop_lost", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Stop Lost"	
		}
		else if(  length( grep( "nonsense", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Nonsense"	
		}
		else if(  length( grep( "missense", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Missense"	
		}
		else if(  length( grep( "ess_splice", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Ess Splice"	
		}
		else if(  length( grep( "splice_acceptor", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Ess Splice"	
		}
		else if(  length( grep( "splice_region", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Splice Region"	
		}
		else if(  length( grep( "nc_variant", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Noncoding Variant"	
		}
		else if(  length( grep( "cds_disrupted", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "CDS Disrupted"	
		}
		else if(  length( grep( "complex_sub", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Complex"	
		}
		else if(  length( grep( "inframe", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Inframe"	
		}
		else if(  length( grep( "nc_transcript_variant", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Noncoding Transcript Variant"	
		}
		else if(  length( grep( "Splice Acceptor", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Ess Splice"	
		}
		else if(  length( grep( "splice_donor", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Splice Donor"	
		}
		else if(  length( grep( "frameshift", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Frameshift"	
		}
		else if(  length( grep( "inframe_deletion", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Inframe Deletion"	
		}
		else if(  length( grep( "synonymous", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Silent"	
		}
		else if(  length( grep( "silent", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Silent"	
		}
		else if(  length( grep( "stop_retained", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Silent"	
		}
		else if(  length( grep( "Del", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Deletion"
		}	
		else if(  length( grep( "Ins", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Insertion"
		}
		else if(  length( grep( "intergenic", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Intergenic"
		}
		else if(  length( grep( "intron_variant", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Intronic"
		}
		else if(  length( grep( "stop_gained_NMD_transcript_variant", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Stop_Gained_NMD"	
		}
		else if(  length( grep( "missense_variant_NMD_transcript_variant", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Missense_NMD"	
		}
		else if(  length( grep( "NMD", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "NMD"	
		}
		else if(  length( grep( "intronic", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Intronic"	
		}
		else if(  length( grep( "UTR", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "UTR"	
		}
		else if(  length( grep( "upstream", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Upstream"	
		}
		else if(  length( grep( "downstream", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Downstream"	
		}
		else if(  length( grep( "non_coding", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Noncoding"	
		}
		else if(  length( grep( "initiator_codon_variant", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Initiator Codon"	
		}
		else if(  length( grep( "initiator_codon_change", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Initiator Codon"	
		}
		else if(  length( grep( "mature_miRNA_variant", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Mature miRNA"	
		}
		else if(  length( grep( "mature miRNA", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Mature miRNA"	
		}
		else if(  length( grep( "splice_region_variant", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Splice Region"	
		}
		else if(  length( grep( "start_lost", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Start Lost"	
		}
		else if(  length( grep( "stop_lost", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Stop Lost"	
		}
		else if(  length( grep( "splice_site_variant", variant_conseq[i] ) ) > 0 )
		{
			variant_conseq_remapped[i] = "Splice Site Variant"	
		}
		else{
			variant_conseq_remapped[i] = variant_conseq[i]
		}
	}

	return(variant_conseq_remapped);
}


## Drop rows and columns from from matrix ##
##
##	Input ::
##	1. Data Matrix
##	2. Row Names to drop
##	3. Col Names to drop

drop_row_cols = function( Matrix, rows, cols )
{
	rows_to_drop = -which( rownames(Matrix) %in% rows )
	if( length( rows_to_drop ) > 0 )
	{
		Matrix_Cleaned = Matrix[ rows_to_drop , ]
	}else{
		Matrix_Cleaned = Matrix
	}
	cols_to_drop = -which( colnames( Matrix ) %in% cols )
	if( length( cols_to_drop ) > 0 )
	{
		Matrix_Cleaned = Matrix_Cleaned[ , cols_to_drop ]
	}else{
		Matrix_Cleaned = Matrix_Cleaned
	}

	return( Matrix_Cleaned );
}

seq_complement <- function( Seq = NULL )
{
  mapping = c("A-T","G-C");
  matched = mapping [ grep( Seq, mapping ) ];
  matched = gsub("-","",matched,fixed=TRUE);
  matched = gsub(Seq,"",matched,fixed=TRUE);
  return(matched);
}

## ---- Sequence Reverse Compliment --- ##
reverse_complement <- function( Seq = NULL )
{
	Seq_Vector = unlist(strsplit(Seq,split=""))
	Seq_Vector_Rev = rev( Seq_Vector );
	Seq_Vector_Rev_Compl = sapply(Seq_Vector_Rev,seq_complement)
	Seq_Rev_Compl = paste( Seq_Vector_Rev_Compl, collapse = "" );
	return(Seq_Rev_Compl)	
}

## -------- ##
#Function :: Maps 12 alterations to 6 
##-------- ##
alteration_remapper = function( alterations = NULL )
{
	
	remapped_alterations = c();
	for( i in 1:length(alterations) )
	{
		if( alterations[i] == "G->A" || alterations[i] == "C->T" )
		{
			remapped_alterations = append( remapped_alterations , "C->T" );	
		}
		else if( alterations[i] == "G->T" || alterations[i] == "C->A" )
		{
			remapped_alterations = append( remapped_alterations , "C->A" );	
		}
		else if( alterations[i] == "G->C" || alterations[i] == "C->G" )
		{
			remapped_alterations = append( remapped_alterations , "C->G" );	
		}
		else if( alterations[i] == "A->G" || alterations[i] == "T->C" )
		{
			remapped_alterations = append( remapped_alterations , "T->C" );	
		}
		else if( alterations[i] == "A->C" || alterations[i] == "T->G" )
		{
			remapped_alterations = append( remapped_alterations , "T->G" );	
		}
		else if( alterations[i] == "A->T" || alterations[i] == "T->A" )
		{
			remapped_alterations = append( remapped_alterations , "T->A" );	
		}
		else{
			remapped_alterations = append( remapped_alterations , "Unrecognized Alteration" );
		}
	}
	
	return( remapped_alterations );
	
}    

alteration_remapper1 = function( alterations = NULL )
{
	
	remapped_alterations = c();
	for( i in 1:length(alterations) )
	{
		if( alterations[i] == "G>A" || alterations[i] == "C>T" )
		{
			remapped_alterations = append( remapped_alterations , "C>T" );	
		}
		else if( alterations[i] == "G>T" || alterations[i] == "C>A" )
		{
			remapped_alterations = append( remapped_alterations , "C>A" );	
		}
		else if( alterations[i] == "G>C" || alterations[i] == "C>G" )
		{
			remapped_alterations = append( remapped_alterations , "C>G" );	
		}
		else if( alterations[i] == "A>G" || alterations[i] == "T>C" )
		{
			remapped_alterations = append( remapped_alterations , "T>C" );	
		}
		else if( alterations[i] == "A>C" || alterations[i] == "T>G" )
		{
			remapped_alterations = append( remapped_alterations , "T>G" );	
		}
		else if( alterations[i] == "A>T" || alterations[i] == "T>A" )
		{
			remapped_alterations = append( remapped_alterations , "T>A" );	
		}
		else{
			remapped_alterations = append( remapped_alterations , "Unrecognized Alteration" );
		}
	}
	
	return( remapped_alterations );
	
}    


## Mapping Base Chagnge to Tri-Nucleotide Context ##

base_change_mapping = function(x,base_changes)
		{
			remapped_base_change = c();
			for( i in 1:length(x) )
			{
				idx = grep( x[i] , base_changes );
				print(i)	
				print(idx)
				remapped_base_change = append( remapped_base_change, base_changes[idx] )
			}		
			return(remapped_base_change)
		}
