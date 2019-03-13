## RScript : PCA.R
## PCA for Dimensionaly reductiona and Principle Component Identification ##
##
## Parameters :: 
##	1. Data Matrix File 
##  2. Plot Directory 
##	3. Label = TRUE | FALSE
##	4. Label Column = 36
##  5. Data Start Column = 8
##  6. Data End Column = 35
##	7. Adjust type = 0 | 1 | 2
##	8. PCA type = sample [row] | feature [column] | both
##	9. skip_rows = 0
## 10. skip_columns = 1
##
##	Example Run : /software/R-2.15.2/bin/Rscript /nfs/users/nfs_r/rm8/My_Code/R/Data_Analysis_Functions/PCA.R /lustre/scratch110/sanger/rm8/experiments/Mutation_Pattern/Data_Sets/Ludmils_Data/Mutation_Data/Annotation/Phase_1_ICGC_WGS_NCSNV_W1_R3/Files/ICGC_W1R3_D_W100R2_P_Tr_Feature_Data_DataFrame.txt /lustre/scratch110/sanger/rm8/experiments/Mutation_Pattern/Data_Sets/Ludmils_Data/Mutation_Data/Annotation/Phase_1_ICGC_WGS_NCSNV_W1_R3/Plots TRUE 182 1 181 0 0 1
##

#Data_File = "Test.txt";
#out_dir = "/Users/rm8/Sanger/NGS_Analysis/experiments/Human/Marco_Cell_Line_PCA/Plots/";
#label = FALSE;
#label_column = 182;
#data_start_column = 1;
#data_end_column = 181;
#adjust_type = 0;
#pca_type = "feature";
#skip_rows = 0;
#skip_columns = 1;


args = commandArgs(TRUE);
Data_File = args[1];
out_dir = args[2];
label = args[3];
label_column = as.numeric(args[4]);
data_start_column = as.numeric(args[5]);
data_end_column = as.numeric(args[6]);
adjust_type = args[7];
pca_type = args[8];
skip_rows = args[9];
skip_columns = args[10];


library(gridExtra)
library(ggplot2)
library(data.table)

print(paste("File : ", Data_File ))
print(paste("Out Dir : ", out_dir ))
print(paste("Label : ", label ))
print(paste("Label Column : ", label_column ))
print(paste("Data Start Column : ", data_start_column ))
print(paste("Data End Column : ", data_end_column ))
print(paste("Adjust Type : ", adjust_type ))
print(paste("PCA Type : ", pca_type ))
print(paste("Skip rows : ", skip_rows ))
print(paste("Skip columns : ", skip_columns ))


source("/nfs/users/nfs_r/rm8/My_Code/R/Data_Processing_Functions/R_Tool_Box.R");
#source("/Users/rm8/Sanger/My_Code/R/Data_Processing_Functions/R_Tool_Box.R")

Data = as.data.frame( fread( Data_File , sep = "\t", header = TRUE, stringsAsFactors = F , skip = skip_rows, check.names=FALSE ) ); 
Data = Data[ , (as.numeric(skip_columns) + 1) : dim(Data)[2] ]

if( label == TRUE )
{
	Data_Label = Data[,label_column]
	color1 = designToColor( Data_Label );
	#pch1 = designToPCH( Data_Label );
	Data = Data[ , data_start_column:data_end_column ]
}

## Here we can take two routes . Either add a noise for columns completely empty
## Or we can remove them from the analysis 


Adjusted_Data = as.matrix(Data)
mode(Adjusted_Data) = "numeric"

if( adjust_type == 0 ) {
	print("No Adjustment done to the data ")
} else if( adjust_type == 1 ) {
	print( "Small Noise will be added to the data" )
	for( i in 1:dim(Adjusted_Data)[2] )
	{
		if( sum ( Adjusted_Data[, i] ) == 0  )
		{
			print (i);
			Adjusted_Data[1, i] = 0.01 ;
		}
	}

	Adjusted_Data = apply ( Adjusted_Data  , 2 , function(X) { ( X - mean(X) ) / sd(X) } );  
} else {
	print( "Rows with no values will removed " )
	
	Adjusted_Data = Adjusted_Data[ rowSums( Adjusted_Data, na.rm = TRUE ) != 0 , ];
	Adjusted_Data = Adjusted_Data[ !is.na( Adjusted_Data[,1] ) ,]
}

print("Data Pre-Processing Done");

# Uncomment this if you want to run it for random 100 features 
# Adjusted_Data = Adjusted_Data[, sample( 1:dim(Adjusted_Data)[2], 100, replace=F )]

if( pca_type == "feature" || pca_type == "both" )
{
	## Feature PCA ##	
	Feature_COV_Matrix =  t( Adjusted_Data ) %*% Adjusted_Data ;  
	Feature_COV_Matrix_EVV = eigen ( Feature_COV_Matrix , symmetric = T );

	Feature_COV_Matrix_EVV_PCS = Feature_COV_Matrix_EVV$vectors[ ,1:10 ];
	Feature_COV_Matrix_EVV_EVALS = Feature_COV_Matrix_EVV$values[ 1:10 ] ;  

	write.table(Feature_COV_Matrix_EVV_PCS,"Feature_COV_Matrix_EVV_PCS.txt",sep="\t",quote=F);

	## PCA  Plot##
	pdf( paste( out_dir, "/Feature_PCA.pdf", sep = "") ,width=18,height=18);
	plot(Feature_COV_Matrix_EVV_PCS[,1],Feature_COV_Matrix_EVV_PCS[,2] , main = "Feaure PCA : PC1 vs PC2" );
	text( Feature_COV_Matrix_EVV_PCS[,1],Feature_COV_Matrix_EVV_PCS[,2] , labels = colnames(Adjusted_Data) );
	dev.off();	

} else if( pca_type == "sample" || pca_type == "both" ) {

	## Mutation PCA ##
	Mutation_COV_Matrix =  Adjusted_Data %*% t(Adjusted_Data) ;  
	Mutation_COV_Matrix_EVV = eigen ( Mutation_COV_Matrix , symmetric = T );

	Mutation_COV_Matrix_EVV_PCS = Mutation_COV_Matrix_EVV$vectors[ ,1:10 ];
	Mutation_COV_Matrix_EVV_EVALS = Mutation_COV_Matrix_EVV$values[ 1:10 ] ;  

	write.table(Mutation_COV_Matrix_EVV_PCS,"Mutation_COV_Matrix_EVV_PCS.txt",sep="\t",quote=F);

	mut = data.frame( PCS=Mutation_COV_Matrix_EVV_PCS, label=Data_Label )
	pdf( paste( out_dir, "/Mutation_PCA_NoJitter.pdf", sep = ""), width=20,height=20);
	p1 = ggplot(mut, aes(x=PCS.1, y=PCS.2 ) ) + geom_point(aes(colour=factor(label)));
	p2 = ggplot(mut, aes(x=PCS.1, y=PCS.3 ) ) + geom_point(aes(colour=factor(label)));
	p3 = ggplot(mut, aes(x=PCS.2, y=PCS.3 ) ) + geom_point(aes(colour=factor(label)));
	p4 = ggplot(mut, aes(x=PCS.2, y=PCS.4 ) ) + geom_point(aes(colour=factor(label)));
	grid.arrange(p1, p2, p3, p4, ncol=2,  main="Feature PCA")
	dev.off();	

} else {
	print("Please check PCA Type")
}

	
## --- Save R Data Object ---- ##
#save.image( paste( out_dir, "/PCA.RData", sep = "" ) );
