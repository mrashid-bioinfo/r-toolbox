## --------------------------------------------------------------------- ##
##					
##				Sub sample using categorical stratification

require(data.table)
require("survey")
data("mtcars")


## ------------------------------------------------- ##
##		Function to stratified sample selection
##
## function to take a random proportional stratified sample of size n
rpss <- function(stratum, n) 
{
    props <- table(stratum)/length(stratum)
    nstrat <- as.vector(round(n*props))
    nstrat[nstrat==0] <- 1
    names(nstrat) <- names(props)
    stratsample(stratum, nstrat)
}


## summary of old datasource ##
# phenodata_summary = function(df)
# {
# 	gender_df = as.data.frame( df %>% group_by( gender, label ) %>% summarise( count = n() ) )
# 	age_df = as.data.frame( df %>% group_by( Age, label ) %>% summarise( count = n() ) )
# 	vendor_df = as.data.frame( df %>% group_by( Vendor, label ) %>% summarise( count = n() ) )
# 	all_summary = rbind( gender_df , age_df, vendor_df )
# 	return( all_summary )
# }
	
## Example Train test data creation ##
create_train_test = function( df, number_of_samples = 20, interaction_vars = c("gear", "carb") )
{
	selrows_train <- rpss( stratum = interaction(df[,interaction_vars], drop=TRUE), n = number_of_samples )

	train_df = df[selrows_train, ]
	write.table( train_df, "train_set.tsv", sep = "\t", quote = F, row.names = F );
	test_df = df[ -selrows_train, ];
	write.table( test_df, "test_set.tsv", sep = "\t", quote = F, row.names = F );
}


args <- commandArgs(trailingOnly = TRUE)
print(args)
print(length(args));

if( length(args) < 3 )
{
    cat("Not enough paramters passed. You need to pass the following paramters \n\n");
    cat(paste(" List of paramters : ",  "\n"));
    cat("\n================================")
    cat("\n1. Dataframe txt file")
    cat("\n2. number of items / samples to select")
    cat("\n3. interaction variables")
    command = paste0(" \n Example Usage  ::  \n")
    cat( command );
    quit(save = "no", status = 0, runLast = TRUE);
}else{

    print( "Found necessary inputs " )
    df_file = args[1];
    no_samples = as.numeric(args[2]);
    interaction_vars = unlist( strsplit( args[3], split = ",", fixed = TRUE ) )
    
    df = read.delim( df_file, sep = "\t", stringsAsFactors = F )
    print(head(df))
    create_train_test(df = df, number_of_samples = no_samples, interaction_vars = interaction_vars )

}

