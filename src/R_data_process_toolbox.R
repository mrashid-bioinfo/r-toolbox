
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#	Given multiple factors the function adds a numeric sequence 
#	after the character values
#	
# 	Example : 
#		Input : CRC, HV, CRC, HV
#		Output : CRC_1, HV_1, CRC_1, HV_1
#
#
#

char_to_num_padded = function(char = NULL)
{
	char_to_num_df = data.frame( char = as.character(char), num = 0 )
	uniq_stage = unique( as.character( char ) )
	
	for( i in 1:length( uniq_stage ) )
	{
		temp_df = char_to_num_df[ char_to_num_df$char ==  uniq_stage[i],  ] 
		char_to_num_df[ char_to_num_df$char ==  uniq_stage[i], "num" ] =  paste( uniq_stage[i], 1:dim(temp_df)[1], sep = "_" )
	}	  

	return( char_to_num_df$num )
}
