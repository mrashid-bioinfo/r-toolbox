## ==================================================================================== ##
##						Data Analysis Toolbox 
##
##


## ==================================================== ##
##
##	Function :: tri_nucleotide_context_fill
##
##	Fill up trinucleotide context
##
##	Paramters 
##	1. Dataframe 
##

tri_nucleotide_context_fill = function( df = NULL , chr = "Chromosome", start = "Start", end = "End", ref = "Reference", alt = "Variant", mutSamples = "mutSamples", species = 'human' )
{
	require(SomaticSignatures)
	require("deconstructSigs")
	require(dplyr)
	require(VariantAnnotation)
	
    df_VR = VRanges( seqnames = df[,chr], ranges = IRanges(start = df[, end], end = df[, end]), ref = df[,ref], alt = df[,alt] )
    
    if( species == "human" )
    {
    	require(BSgenome.Hsapiens.1000genomes.hs37d5)
    	print("Human")	
    	df_VR_context = as.data.frame( mutationContext(df_VR, BSgenome.Hsapiens.1000genomes.hs37d5) )
    }
    else if( species == "mouse" ) {
    	print("Mouse")
    	require(BSgenome.Mmusculus.UCSC.mm10)
    	df_VR_context = as.data.frame( mutationContext(df_VR, BSgenome.Mmusculus.UCSC.mm10) )	
    }
    else{
    	print( "Unknown species" )
    }
    
    print("Creating modified VR object")
    df_VR_context_mod = df_VR_context %>% separate(context, c("Preceeding","Trailing") )
    df_VR_context_mod1 = df_VR_context_mod %>% separate(alteration, c("Ref","Alt"), sep = 1 )
    df_VR_context_mod1$tri_context = paste( df_VR_context_mod1$Preceeding, "[", df_VR_context_mod1$Ref, ">" , df_VR_context_mod1$Alt, "]", df_VR_context_mod1$Trailing, sep = "" )

    print("Creating Data Frame from VR")
    df_VR_context_df = data.frame( Chromosome = df[,chr], Start = df[,start], End = df[,end], Reference = df[,ref], Variant = df[,alt], tri_context = df_VR_context_mod1$tri_context,mutSamples = df[,mutSamples] )

    Spir_triContext_perSample = table( df_VR_context_df$mutSamples, df_VR_context_df$tri_context )
    return(Spir_triContext_perSample)
}


## ==================================================== ##
##
##	Function 	:: SomSig_Analysis
##
##	Description :: Assess COSMIC signature contribution using deConstructSig
##
##	Paramters 
##	1. Dataframe 
##


SomSig_Analysis = function( df = NULL, chr = "Chromosome", start = "Start", end = "End", ref = "Reference", alt = "Variant", mutSamples = "mutSamples", species = 'human', norm_method = 'exome' )
{

	require(SomaticSignatures)
	require("deconstructSigs")
	require(dplyr)
	require(VariantAnnotation)
	

	VR = VRanges( seqnames = df[,chr], ranges = IRanges(start = df[, end], end = df[, end]), ref = df[,ref], alt = df[,alt] , sampleNames = df[,mutSamples] )
	
	if( species == 'human' )
	{	
		require(BSgenome.Hsapiens.1000genomes.hs37d5)
		print("Human")	
		VR_context = mutationContext(VR, BSgenome.Hsapiens.1000genomes.hs37d5)	
	}
	else if( species == "mouse" ) {
		library(BSgenome.Mmusculus.UCSC.mm10)
		print("Mouse")	
		VR_context = mutationContext(VR, BSgenome.Mmusculus.UCSC.mm10)		
	}
	else{
		print( "Unknown species" )
	}
	
	sca_mm = motifMatrix(VR_context, normalize = F)
	mm_sca=as.data.frame(t(sca_mm))
	colnames(mm_sca)=colnames(signatures.cosmic)

	uniq_samples = unique(df[,mutSamples])
	s=c()
	for (i in 1:length(uniq_samples)) 
	{
		test1 = whichSignatures(tumor.ref = mm_sca,signatures.ref =signatures.cosmic,sample.id =uniq_samples[i] , contexts.needed = TRUE, tri.counts.method = norm_method )
		weights1<- data.frame(test1[["weights"]])
		unknown<- test1[["unknown"]]
		weights1$unknown <- unknown
		s=rbind(s,weights1)

	}

  	return(s)
}




## Mean, Minimum and Maximum of data matrix ##
##	
##	Computes mean and standard deviation of a matrix [ by column or row ]
##		
##	Input :: 
##		1. matrix
##		2. direction [ 1 for row; 2 for column ]
#
##	output :: matrix
##
min_max_mean = function( mat = null, direction = 1 )
{
	if( direction == 1 )
	{
		avg = rowMeans( mat , na.rm = TRUE )
		minx = apply( mat , 1, function(X){ min(X, na.rm = T ) } )
		maxx = apply( mat , 1, function(X){ max(X, na.rm = T ) } )
	}
	else if( direction == 2 )
	{
		avg = colMeans( mat , na.rm = TRUE )
		minx = apply( mat , 2, function(X){ min(X, na.rm = T ) } )
		maxx = apply( mat , 2, function(X){ max(X, na.rm = T ) } )
	}
	else{
		print( "Invalid direction : Only 1[row] and 2 [column] is allowed" )
	}

	avg[ is.na(avg) ] = 2
	minx[ is.infinite(minx) ] = 2
	maxx[ is.infinite(maxx) ] = 2
	df = data.frame( avg, minx, maxx );
	return(df)
}



## Mean, SD, Lower and Upper ##
##	
##	Computes mean and standard deviation of a matrix [ by column or row ]
##		
##	Input :: 
##		1. matrix
##		2. direction [ 1 for row; 2 for column ]
#
##	output :: dataframe
##

mean_sd = function( mat = null, direction = 1 )
{
	if( direction == 1 )
	{
	avg = rowMeans( mat , na.rm = TRUE )
	sd = apply( mat , 1, function(X){ sd(X, na.rm = T ) } )
	}
	else if( direction == 2 )
	{
	avg = colMeans( mat , na.rm = TRUE )
	sd = apply( mat , 2, function(X){ sd(X, na.rm = T ) } )
	}
	else{
	print( "Invalid direction : Only 1[row] and 2 [column] is allowed" )
	}
	sd[ is.na(sd) ] = 0.05
	lower = avg - sd
	upper = avg + sd
	df = data.frame( avg, sd, lower, upper );
	return(df)
}



## Bed region overlap ##

find_overlap_region = function( df )
{
	## !!!! §§§ Still incomplete ##	
	for( i in 1:dim(cnv_sample_df_mod)[1] )	
	{
		if( i == 1 )
		{
			old = cnv_sample_df_mod[i,]
		}
		else
		{
			curr = cnv_sample_df_mod[i,]

			if( old[,"chromosome"] == curr[,"chromosome"] )
			{
				if( ( as.numeric( curr[ , "start.pos" ] ) <= as.numeric( old[ , "end.pos" ] ) )  && ( as.numeric( curr[ , "end.pos" ] ) >= as.numeric( old[ , "end.pos" ] ) ) )
				{
					old[,"end.pos"] = curr[,"end.pos"]

				}
			}
		}
	}
}


## Chromosome Rename ##
chr_rename = function( chr )
{
	a = gsub( "X", "23", chr )
	a = gsub( "Y", "24", a );
	return( a )
}


## Sequennza CNV Classification ##
sequenza_cnv_conversion = function( df )
{
	cnv_class = apply( df, 1, function(X)
	{
		if( X[2] < 100 )
		{
			return("Too Small")
		}
		# Focal amps : 
		else if( X[1] >= 5 &&  X[2] <= 1000000 )
		{
			return("Focal Amplification")
		}	
		# Large amps : 
		else if( X[1] >= 5 &&  X[2] > 1000000 )
		{
			return("Large Amplification")
		}
		# Focal Gain : 
		else if( X[1] >= 3 &&  X[2] > 100 && X[2] <= 1000000 )
		{
			return("Focal Gain")
		}
		# Small Gain : 
		else if( X[1] >= 3 &&  X[2] > 1000000 )
		{
			return("Large Gain")
		}
		# Focal Deletion
		else if( X[1] == 0 &&  X[2] > 100 && X[2] <= 1000000 )
		{
			return("Focal Deletion")
		}
		# Deletion
		else if( X[1] == 0 &&  X[2] > 1000000 )
		{
			return("Deletion")
		}
		# Loss
		else if( X[1] == 1 &&  X[2] > 100 )
		{
			return("Loss")
		}
	} )

	return( cnv_class );
}