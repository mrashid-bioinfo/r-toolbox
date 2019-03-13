###################################
#                                 #
####    R Library    ##############
#                                 #
###################################

## --------------------------------------------------------------- ##
##	Function :: Collapse rows of a data frame to one row
##
##		Input Paramter :: 
##		1. df 				= data frame
##

collapse_rows = function( df )
{
	apply( df, 2, function( X )
	{ 
		uniq_x = unique(X);
		uniq_x = uniq_x[ uniq_x!= 0 ];

		if( length( uniq_x ) == 0 )
		{
			x_collapse = "0"
		}
		else{
			x_collapse = paste( uniq_x, sep = "," );
		}

		return( x_collapse )
	} )
}

## --------------------------------------------------------------- ##
##	Function :: Summarise dataframe values to specified intervals
##
##		Use  :: Used to summarise data over an specified interval
##
##		Input Paramter :: 
##		1. df 				= data frame
##		2. interval_col		= Column to summrise over 
##		3. summarise_col	= Column value to summrise  
##		4. vector  			= specified interval
##

summarise_col_over_interval = function( df = NULL, interval_col = NULL, summarise_col = NULL, vector = c( seq(0,100,10), 10000000 ) )
{
	df[,interval_col] = as.numeric(df[,interval_col])
	df[,summarise_col] = as.numeric(df[,summarise_col])

	summarised_df = c()
	for( i in 1:( length(vector)-1 ) )	
	{
		temp = sum( df[ df[,interval_col] > vector[i] & df[,interval_col] <= vector[i+1], summarise_col ] );
		summarised_df = rbind( summarised_df, c( paste( vector[i], vector[i+1] , sep = "-" ) , temp ) );
		print(i)
		print(summarised_df)
	}

	return(summarised_df)
}


## ---------------------------------------------------------- ##
##	Function :: To find end of a genomic chunk
##
##		Use  :: Used in Binomial test wrapper
##
##		Input Paramter :: 
##		1. temp 				= data frame
##		2. end_row  			= index of the last row
##		3. largest_window_size 	= 50
##

find_end = function( temp = NULL, end_row = NULL, largest_window_size = 50 )
{
	counter = 0	
	while( TRUE )
	{
		current_pos = temp[end_row,3]	
		next_pos = temp[end_row+1,3]

		## Make sure you do not split position closer than 50 bp [or whatever the largest window size is] 
		## And this is fixed at largest window size. This do not change with window size . 

		if( ( next_pos - current_pos ) > largest_window_size )
		{
			return(end_row)
		}
		else{
			end_row = end_row + 1	
		}
		counter = counter + 1;
	}
	return(end_row)
}




## CGP mutation data format to Somatic Signature Format  ##
## 

CGP_to_SomaticSig = function( df = NULL, cols = c("Chromosome", "End", "Reference", "Variant"), tag = "PD" )
{
	require( dplyr )
	require( tidyr )

	df$mut_samples = get_mutated_samples( df, tag = tag )
	df_sig_E = as.data.frame( df %>% mutate( mut_samples = strsplit( as.character( mut_samples ), split = "," ) ) %>% unnest( mut_samples ) )
	df_sig_E1 = as.data.frame( df_sig_E %>% mutate( Variant = strsplit( as.character( Variant ), split = "," ) ) %>% unnest( Variant ) )

	df_sig = data.frame ( Chromosome = df_sig_E1[, cols[1]], start = df_sig_E1[, cols[2]], end = df_sig_E1[, cols[2]], Reference = df_sig_E1[, cols[3]], Variant = df_sig_E1[, cols[4]],  mut_samples = df_sig_E1$mut_samples )

	return(df_sig_E)
}


## CGP mutation data format to EMu Format  ##
## 

CGP_to_Emu = function( df = NULL, cols = c("Chromosome", "End", "Reference", "Variant"), tag = "PD" )
{
	df$mut_samples = get_mutated_samples( df, tag = tag )
	df_sig_E = as.data.frame( df %>% mutate( mut_samples = strsplit( as.character( mut_samples ), split = "," ) ) %>% unnest( mut_samples ) )
	df_sig_E1 = as.data.frame( df_sig_E %>% mutate( Variant = strsplit( as.character( Variant ), split = "," ) ) %>% unnest( Variant ) )

	df_sig = data.frame ( mut_samples = df_sig_E1$mut_samples, Chromosome = df_sig_E1[ , cols[1] ], Locus = df_sig_E1[ , cols[2] ], Nucleotide_Change =  paste( df_sig_E1[ , cols[3] ], ">", df_sig_E1[ , cols[4] ], sep = "" ) )

	if( tag == "PD" )
	{
		df_sig$Chromosome = gsub( "X", "23", df_sig$Chromosome )
		df_sig$Chromosome = gsub( "Y", "24", df_sig$Chromosome )
	}else
	{
		df_sig$Chromosome = gsub( "X", "20", df_sig$Chromosome )
		df_sig$Chromosome = gsub( "Y", "21", df_sig$Chromosome )
	}

	return(df_sig)
}


## CGP mutation data format to EMu Format  ##
## 

CGP_to_dNdS = function( df = NULL, cols = c("Chromosome", "End", "Reference", "Variant", "Gene_Name"), tag = "PD" )
{
	df$mut_samples = get_mutated_samples( df, tag = tag )
	df_sig_E = as.data.frame( df %>% mutate( mut_samples = strsplit( as.character( mut_samples ), split = "," ) ) %>% unnest( mut_samples ) )
	df_sig_E1 = as.data.frame( df_sig_E %>% mutate( Variant = strsplit( as.character( Variant ), split = "," ) ) %>% unnest( Variant ) )

	if( length( cols ) == 4 )
	{
		df_dnds = data.frame ( sampleID = df_sig_E1$mut_samples, chr = df_sig_E1[ , cols[1] ], pos = df_sig_E1[ , cols[2] ], ref =  df_sig_E1[ , cols[3] ], mut =  df_sig_E1[ , cols[4] ] )	
	}
	else{
		df_dnds = data.frame ( sampleID = df_sig_E1$mut_samples, chr = df_sig_E1[ , cols[1] ], pos = df_sig_E1[ , cols[2] ], ref =  df_sig_E1[ , cols[3] ], mut =  df_sig_E1[ , cols[4] ], gene = df_sig_E1[ , cols[5] ] )		
	}
	

	return(df_dnds)
}




###################################
##
##	Capitalize the first Letter of the word

simpleCap1 <- function(x) {
  paste( toupper(substring(x,1,1)), substring(x,2,nchar(x)), sep = "")
}

###################################
##
##	Capitalize the first Letter of the word

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}



###################################
## 	Function :: drop_samples
## 	Drop samples from data frame
##
##	Input :: 
##		1. Data Frame [ Alt_Base, foramt field data from CGP VCF in following format :: %PM|%FAZ|%FCZ|%FGZ|%FTZ|%RAZ|%RCZ|%RGZ|%RTZ ]
##		2. all sample ids
##		3. samples to drop
##
##	Output :: 
##		Data frame after dropping sample
##

drop_samples = function( df, sample_ids, sample_ids_to_drop )
{
	dropped_sample_idx = which ( colnames(df) %in% sample_ids_to_drop )
	temp_data = as.matrix(df [ , sample_ids [ ! sample_ids %in% sample_ids_to_drop ] ])
	
	dropped_df = df [ rowSums( temp_data ) > 0 , -dropped_sample_idx ] 

	return( dropped_df )
}


###################################
## 	Function :: mut_read_count
## 	Extract the mutant read count from format field of the 
##
##	Input :: 
##		1. Data Frame [ Alt_Base, foramt field data from CGP VCF in following format :: %PM|%FAZ|%FCZ|%FGZ|%FTZ|%RAZ|%RCZ|%RGZ|%RTZ ]
##		2. Alt base col ID
##		3. Sample start col
##		3. Sample end col
##
##	Output :: 
##		A vector with sum of forward and reverse reads supporting the ALT base 
##


mut_read_count = function( df = NULL, alt_col = 4, sample_start = 11, sample_end = 25 )
{
	
	mut_read_ct = apply( df[,c(alt_col,sample_start:sample_end)], 1 , function(X)
	{
		if( X[1] == "A" )
		{
			mut_rct = sapply( X[2:length(X)], function(Y)
			{ 
				temp = unlist ( strsplit( Y, split = "|", fixed = T ) ); 
				if( temp[1] != "." && temp[5] != "." ) 
				{ 
					sum( as.numeric(temp[1]) + as.numeric(temp[5]) )
				} 
			} )	
		}
		if( X[1] == "C" )
		{
		
			mut_rct = sapply( X[2:length(X)], function(Y)
			{ 
				temp = unlist ( strsplit( Y, split = "|", fixed = T ) ); 
				if( temp[2] != "." && temp[6] != "." ) 
				{ 
					sum( as.numeric(temp[2]) + as.numeric(temp[6]) )
				} 
			} )	
		}
		if( X[1] == "G" )
		{
			mut_rct = sapply( X[2:length(X)], function(Y)
			{ 
				temp = unlist ( strsplit( Y, split = "|", fixed = T ) ); 
				if( temp[3] != "." && temp[7] != "." ) 
				{ 
					sum( as.numeric(temp[3]) + as.numeric(temp[7]) )
				} 
			} )	
		}
		if( X[1] == "T" )
		{
			mut_rct = sapply( X[2:length(X)], function(Y)
			{ 
				temp = unlist ( strsplit( Y, split = "|", fixed = T ) ); 
				if( temp[4] != "." && temp[8] != "." ) 
				{ 
					sum( as.numeric(temp[4]) + as.numeric(temp[8]) )
				} 
			} )	
		}
	
		max( unlist( mut_rct ) );
	} )

	return(mut_read_ct);
}
		

###################################
## 	Function :: merge_by_col
## 	Merge [ Both sum and Mean ] a data frame based on one column id
##
##	Input :: 
##		1. Data Frame
##		2. Column ID
##		..
##		..
##
##	Output :: 
##		List containing two data frames
##

merge_by_col = function( df = NULL, id_cols = 1:23, feature_start = 24, feature_end = 218, Data_Label = 219, segment_label = "Genome_Segment_Label_1K" )
{
	uniq_segs = unique( df$Genome_Segment_Label_1K )
	for( i in 1:length( uniq_segs ) )
	{
		temp = df[ df[,segment_label] == uniq_segs[i], ]

		temp_ids = apply( temp[ , id_cols], 1, function(X){
			paste( X, collapse = "|" )
		} )

		merged_ids = paste( temp_ids, collapse = "%" );
		merged_label = sum( temp[,Data_Label] );
		sum_data = colSums( temp[,feature_start:feature_end] );
		mean_data = colMeans( temp[,feature_start:feature_end] );

		temp_sum_df = c( Merged_ID = merged_ids, Merged_Label = merged_label, sum_data )
		temp_mean_df = c( Merged_ID = merged_ids, Merged_Label = merged_label, mean_data )

		if( i == 1 )
		{
			sum_df = temp_sum_df
			mean_df = temp_mean_df	
		}
		else
		{
			sum_df = rbind( sum_df, temp_sum_df );
			mean_df = rbind( mean_df, temp_mean_df );
		}
	}

	return( list( sum_df, mean_df ) );
}
		

## Sigmoid function to scale vector 
sigmoid_np = function(x)
{
    	return ( 1 / ( 1 + exp(-x) ) );	
}

## Function to subsample from two vectors ##
vector_pair_subset = function( x, y, size = 100 )
	{
		size = size - 2
		x_1	 = x[1]
		y_1 = y[1]

		idx = sort( sample( 2:( length( x ) - 1) , size , replace = F ) ) ;
		x_1 = append( x_1 , x[ idx ] )
		y_1 = append( y_1 , y[ idx ] )
		
		x_1 = append( x_1, x[ length(x) ] )
		y_1 = append( y_1, y[ length(y) ] )

		return( list ( x_1 , y_1 ) );
	}


## Function to return static 96 mutation context ##

mutation_96_context = function( sep = "bracket" )
{
	bases = c('A','C','G','T');
	core_base = c( "C>A","C>G","C>T", "T>A", "T>C", "T>G" )
	base_96 = c();
	for( i in 1:length(core_base) )
	{
		for(j in 1:length(bases))
		{
			for(k in 1:length(bases))
			{
				if(sep == "bracket")
				{
					base_96 = append( base_96, paste( bases[j], "[", core_base[i], "]", bases[k] , sep = "" ) );	
				}
				else{
					base_96 = append( base_96, paste( bases[j], " ", core_base[i], " ", bases[k] , sep = "" ) );	
				}
			}
		}
	}
	return (base_96)
}



## Function to split 2 level INFO  ##

INFO_Split = function( INFO, level1_tag = "CSQ", level2_index = 4 )
{
        unlist ( lapply( as.list ( INFO ) , function(X)
        {
                split1 = unlist ( strsplit(X, split = ";" ) )
                a = split1 [ grep( level1_tag, split1 ) ]
				b = unlist ( strsplit(a, split = "|", fixed = T ) ) [level2_index]
                return(b)
        } ) )
}

## Function to rowmax from Data Frame / matrix ##
rowMax = function( x, na.rm = T  )
{
    apply(x,1,function(X)
    {
        max(X, na.rm = na.rm)
    })
}

## Function to rowmax from Data Frame / matrix ##
colMax = function( x, na.rm = T  )
{
    apply(x,2,function(X)
    {
        max(X, na.rm = na.rm)
    })
}


## Function to rowmins from Data Frame / matrix ##
rowMins = function( x, na.rm = T )
{
    apply(x,1,function(X)
    {
        min(X, na.rm = na.rm )
    })
}

## Function to rowSD from Data Frame / matrix ##
rowSD = function( x, na.rm = T )
{
	print("Warning. In case you supplied single value, \n rowSD function will return 0 as SD. ")
    apply(x,1,function(X)
    {
        sd_val = sd(X, na.rm = na.rm )
		if( is.na(sd_val))
        {
        	return(0)	
        }
        else{
        	return(sd_val)	
        }
    })
}

## Function to row Meadia from Data Frame / matrix ##
rowMedian = function( x, na.rm = T  )
{
    apply(x,1,function(X)
    {
        median(X, na.rm = na.rm)
    })
}

## Function to col Meadia from Data Frame / matrix ##
colMedian = function( x, na.rm = T  )
{
    apply(x,2,function(X)
    {
        median(X, na.rm = na.rm)
    })
}


## Function to Make ID from Data Frame ##
return_id_stac = function( df )
{
    id = paste( "chr", df[,1], ":" , df[,2], "-", df[,3], sep = "" );
    return(id);
}


return_id = function( df, cols = c(1,2,3,4,5), collapse = "-" )
{
	return( apply( df[,cols], 1, function(X){ X = gsub("\\s+", "", X ) ; return(paste( as.character(X), collapse = collapse) ) } ) )
}


## Function to Make ID from Data Frame ##
# return_id = function( df, count = 5 )
# {
# 	if( count == 5 )	
# 	{
# 		id = paste( df[,1], df[,2], df[,3], df[,4], df[,5], sep = "_" );
# 	}
# 	if( count == 4 )		
# 	{
# 		id = paste( df[,1], df[,2], df[,3], df[,4], sep = "_" );
# 	}
# 	if( count == 3 )
# 	{
# 		id = paste( df[,1], df[,2], df[,3], sep = "_" );	
# 	}
# 	if( count == 2 )
#         {
#                 id = paste( df[,1], df[,2], sep = "_" );
#         }
# 	return(id);
# }


##
## Brings copy number data between [ -2  -- 0 --- 10 Range]
clip_cnv = function(dat=NULL,base_copy=2,bottom=-2,top=10,data_start_col=4,data_end_col=39)
{
    apply( dat[,data_start_col:data_end_col], 2, function(X)
    {
        X = X - base_copy
        X = ifelse( X >= top, top, ifelse( X <= bottom, bottom, X ) )
    }
    )
}

 most_common_value <- function(x)
 {	 
	 return( as(names(which.max(table(x))), mode(x)) ); 
 } 
  
 recurrent_value_count <- function(x)
 {	 
	 a = sort(table(x)) 
 }

## ------- ##
## Function to create test and training data 
## ------- ##
 
split_training_test = function( data_frame = NULL, fraction = 0.5 )
		{
			split_index_size = round ( dim(data_frame)[1] * fraction );
			set.seed(10);
			split_index = sample( dim(data_frame)[1], split_index_size , replace = F );

			tr = data_frame[ split_index, ]
			te = data_frame[ -split_index, ]

			return( list( tr, te ) );
		}


## ----- ##
## list_to_vector : Function to extract a specific entry/pattern from a composite list and return as a vector ## 
## ----- ##

list_to_vector <- function( list_l = NULL , list_entry = NULL , pattern = NULL)
{
	a=c(); 
	
	if( !is.null(list_entry) ){
		for(i in 1:length(list_l)) { a = append(a, list_l[[i]][list_entry]) }	
	}
	else{
		
		if( !is.null(pattern) ){
			for(i in 1:length(list_l)) { a = append(a, list_l[[i]][grep( pattern, list_l[[i]] , fixed = TRUE )] ) }
		}
		else{
			print("None of the paramter is defined . will return an empty vector");
		}
	}
	
	return(a);
}

## ----- ##
## vector_to_matrix : Function to split a vector and return as matrix ## 
## ----- ##
vector_to_matrix <- function( vector = NULL , split = "", mat_col = NULL )
{
	a = unlist(strsplit(vector,split=split));
	dim(a) = c( mat_col, length(vector) );
	a = t(a);
	return(a);
}

## ----- ##
## binMatrix : function to bin values of x (vector) based on it's corresponding value in y
## ----- ##
 
binMatrix <- function( x = NULL , y = NULL, no_of_bins = 10 )
{	
	groups = levels(as.factor(x));
	print(groups);
	bin_matrix = matrix( 0, no_of_bins, length( groups) );
	#print(bin_matrix);
	
	temp = cbind(x,y);
	
	step_breaks = ( max( y ) - min( y ) ) / no_of_bins;
	#print(step_breaks);
	rowname_bin_matrix = c();
	
	for( i in 0:(no_of_bins-1) )
	{
		start = i * step_breaks ;
		end = ( i * step_breaks ) + step_breaks ; 
		
		rowname_bin_matrix = append( rowname_bin_matrix , paste( start, end , sep = " - " ) );
		
		if( i == (no_of_bins-1) )
		{
			end = end + step_breaks;	
		}
		
		subTable = table( temp[ (temp[,2] >= start) & (temp[,2] < end), 1 ] ) ;
		#print(subTable);
		
		for( j in 1:length(groups) )
		{
			if( !is.na( subTable[ groups[j] ] ) ){
				bin_matrix[ (i+1) , j ] = subTable[ groups[j] ];	
			}
			else{
				bin_matrix[ (i+1) , j ] = 0;
			}
			 			
		}	
	}
	
	rownames(bin_matrix) = rowname_bin_matrix;
	colnames(bin_matrix) = groups;
	print(bin_matrix);
	return(bin_matrix);
}


## ---- ## 
## Clip Values to Quantiles ##
## ---- ##

clip_vector_to_quantile <- function( x = NULL , quant = (0:10) * 0.1 , lowest_cutoff = 2, highest_cutoff = 10 )
{
	q_x = quantile(x, quant );
	print(q_x);
	x[which( x <= q_x[ lowest_cutoff ] )] = q_x[lowest_cutoff];
	x[which( x >= q_x[ highest_cutoff ] )] = q_x[highest_cutoff];
	
	return(x);
}

## ------- ##
## Function : Calculate TPR and FPR 
## ------- ##

tpr_fpr_calculation <- function( x = NULL, labels = NULL, p = NULL, n = NULL )
{
	data = cbind( x , labels );
	
	p_n_t = table( labels );
	
	p = p_n_t["1"];
	n = p_n_t["0"];
	
	print(p);
	print(n);
	
	tp_fp_Matrix = matrix( 0, length(x) , 2 );
	
	x_s = sort(x);
	print(x_s);
	for( i in 1:length(x_s) )
	{
		a_t = table( data[ data[,1] <= x_s[i] , 2 ] ) ;
		
		print(a_t);
		
		if( !is.na(a_t["0"]) ) { tp_fp_Matrix[ i, 2 ] = a_t["0"] / n };
		if( !is.na(a_t["1"]) ) { tp_fp_Matrix[ i, 1 ] = a_t["1"] / p };
	}
	
	return( tp_fp_Matrix );
	
}




## ------- ##
## Function : Trim numerics in gene name
## ------- ##

trim_genename = function( genes = NULL)
{
	genes_1 = gsub("[0-9][a-z]","",genes)
	genes_2 = gsub("[0-9]","",genes_1)
	return(genes_2)
}


## ------- ##
## Function : Read delimeted File  
## ------- ##

read_file = function( file_path = NULL, sep = "\t" , header = T, stringsAsFactors = F , row.names = NULL )
{
	file_data = read.delim( file_path , sep = sep, header = header, stringsAsFactors = stringsAsFactors , row.names = row.names );
	#colnames(file_data) = gsub( "X", "", colnames(file_data) );
	return(file_data);
}

## ------- ##
## Function : Mean Normalize a vector  
## ------- ##

normalize_and_scale = function( vector = NULL , bottom = 0 , top = 10  )
{
	
	vector = as.numeric(vector);
	norm_vector = ( vector - mean( vector ) ) / sd( vector ) ;
	min = min(norm_vector);
	norm_scaled_vector = norm_vector - (min - bottom);
	norm_scaled_vector = norm_scaled_vector * ( top / max(norm_scaled_vector) );
	
	return(norm_scaled_vector);
}

## ------- ##
## Function : Rescale vector to a specified range  
## ------- ##

rescale = function(data, bottom = 0, top = 20){
	min1 = min(data);
	max1 = max(data);
	
	scaled_data = ( data - min1 ) *  ( ( top - bottom ) / ( max1- min1 ) ) +  bottom;
	return(scaled_data);
	
}


## Design to color pallet conversion conversion ##

## Inputs : 

# - Expects a vector of design
# - user can also specify colors in order they want
# - Returns a color pallete coded according to factors in design

designToColor = function(design,colorBox= c("red4","chartreuse4","chocolate","cornflowerblue","cornsilk3","gray20","gray50",  "mediumvioletred","orange") , factor_length = 9 )
{
	design <- as.factor(design)
	design_length = length(levels(design));
	if( design_length > factor_length ){    
		colorBox = colors()[ sample(657, design_length,replace=F) ];
	}
	
	#print(design);
	color <- c()
	for(k in 1:length(design)) 
	{ 
		for(l in 1:length(levels(design))) 
		{ 
			if (design[k] == levels(design)[l]) 
				color <- append(color,colorBox[l]); 
		} 
	}
	print(" Return Color ");
	return(color);  
}

## Accepts a vector with values from different class and returns a corresponding R PCH value ##

designToPCH = function(design,pchBox= c(0,1,2,3,4) , factor_length = 9 )
{
	design <- as.factor(design)
	design_length = length(levels(design));
	if( design_length > factor_length ){    
		pchBox = sample(1:18, design_length,replace=F);
	}
	
	#print(design);
	pch <- c()
	for(k in 1:length(design)) 
	{ 
		for(l in 1:length(levels(design))) 
		{ 
			if (design[k] == levels(design)[l]) 
				pch <- append(pch,pchBox[l]); 
		} 
	}
	print(" Return PCH ");
	return(pch);    
}


## List append ##

lappend <- function(lst, obj) {
	lst[[length(lst)+1]] <- obj
	return(lst)
	}

	
## Multiple GGplot plots in one page ##

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
	require(grid)
	
	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)
	
	numPlots = length(plots)
	
	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) {
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
				ncol = cols, nrow = ceiling(numPlots/cols))
	}
	
	if (numPlots==1) {
		print(plots[[1]])
		
	} else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		
		# Make each plot, in the correct location
		for (i in 1:numPlots) {
			# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
			
			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
							layout.pos.col = matchidx$col))
		}
	}
}

## Function to normalize data by mean and standard deviatiion ## 
normalize = function( vector = NULL )
{
	
	print("## Function to normalize data by mean and standard deviation ##");
	
	vector = as.numeric(vector);
	norm_vector = ( vector - mean( vector ) ) / sd( vector ) ;
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
	require(grid)
	
	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)
	
	numPlots = length(plots)
	
	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) {
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
				ncol = cols, nrow = ceiling(numPlots/cols))
	}
	
	if (numPlots==1) {
		print(plots[[1]])
		
	} else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		
		# Make each plot, in the correct location
		for (i in 1:numPlots) {
			# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
			
			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
							layout.pos.col = matchidx$col))
		}
	}
}

## Compares the rowSums for a set of user specified columns and 
## returns the number of rows passed the threshold. 

rowSumCheck = function( column_names, data, value )	
{
	if( length(column_names) > 1 )
	{	
		a = data [ rowSums(data[, column_names ] ) > value , ];
		return( dim(a)[1] / dim(data)[1] );
	}
	else{
		return( dim( data[ data[,column_names] > value, ]  )[1]  / dim(data)[1] );	
	}
}	

## --------- ##
## Count number of non-zero entries
rowCounts = function( mat )
{
	apply( mat, 1, function( X) { length( X[ X > 0 ] ) } )
}

