## ============================================================================ ##
## 			 Function ::  Runs Classifier for multiple negative data size
##		
## ============================================================================ ##

## Function to demonstrate PCA vs tSNE
##
##

pca_vs_tsne = function()
{

	library(readr)
	library(Rtsne)
	# The competition datafiles are in the directory ../input
	# Read competition data files:
	train <- read_csv("../input/train.csv")
	test <- read_csv("../input/test.csv")
	train$label <- as.factor(train$label)

	# shrinking the size for the time limit
	numTrain <- 10000
	set.seed(1)
	rows <- sample(1:nrow(train), numTrain)
	train <- train[rows,]
	# using tsne
	set.seed(1) # for reproducibility
	tsne <- Rtsne(train[,-1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
	# visualizing
	colors = rainbow(length(unique(train$label)))
	names(colors) = unique(train$label)
	plot(tsne$Y, t='n', main="tsne")
	text(tsne$Y, labels=train$label, col=colors[train$label])

	# compare with pca
	pca = princomp(train[,-1])$scores[,1:2]
	plot(pca, t='n', main="pca")
	text(pca, labels=train$label,col=colors[train$label])

	# Generate output files with write_csv(), plot() or ggplot()
	# Any files you write to the current directory get shown as outputs

}


## Function to convert values to nearest quantile ( lower quantile )
##
##  x 		:: numeric value
##	quants 	:: quantile values of vector where x originated from
## 	direction :: greater ( positive values ) or smaller ( negative values )
##

quantile_transform = function( x = NULL, quants = NULL, direction = "greater" )
{
	if( direction == "greater" )
	{
		return( quants [ which ( quants >= x)[1] ] )
	}
	else if( direction == "smaller" )
	{
		return( quants [ tail( which ( quants <= x), 1) ] )	
	}
	else{
		print("Unknown direction")
	}
}

## Function to compute simple AUC
## Function :: simple_auc
##
simple_auc <- function(TPR, FPR)
{
	  # inputs already sorted, best scores first 
	  dFPR <- c(diff(FPR), 0)
	  dTPR <- c(diff(TPR), 0)
	  sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}

##  Function :: auc_range
##  Function to compute simple AUC between range
##
##	Function to compute AUC between a range. Default is between min and max of x

auc_range <- function(x, y, x_start = min(x), x_end = max(x), method = "spline")
{
	require(MESS)
	auc(x,y, from =  x_start, to = x_end, type = method )  
}


## Function :: average_prc_curve 
##	
##	Description : Computes average [average of folds] precision-recall curve using precrec package
##

average_prc_curve = function( df, score_cols = "asym_RFC", lab_col = "bfRec_Label", fold_col = "Fold" )
{
	require("precrec");
	smmdatlabs <- mmdata( nfold_df = df, score_cols = score_cols, lab_col = lab_col, fold_col = fold_col ) 
	smcurveslabs <- evalmod(smmdatlabs, raw_curves = TRUE)
	asym_avg_attr = attr(smcurveslabs$prcs,"avgcurves")

	return(asym_avg_attr)
}

## Function :: average_roc_curve 
##	
##	Description : Computes average [average of folds] receiver operator curve using precrec package
##

average_roc_curve = function( df, score_cols = "asym_RFC", lab_col = "bfRec_Label", fold_col = "Fold" )
{
	require("precrec");
	smmdatlabs <- mmdata( nfold_df = df, score_cols = score_cols, lab_col = lab_col, fold_col = fold_col ) 
	smcurveslabs <- evalmod(smmdatlabs, raw_curves = TRUE)
	asym_avg_attr = attr(smcurveslabs$rocs,"avgcurves")
	return(asym_avg_attr)
}


## Function :: performance_curve 
##	
##	Description : A wrapper function for average_roc and average_prc curve
##

performance_curve = function( df = NULL,  clf_list = clf_list, score_cols = NULL, metric = "prc", lab_col = "bfRec_Label", fold_col = "Fold" )
{
	for( k in 1:length( clf_list ) )
	{
		if( metric == "roc" )	
		{
			attr = average_roc_curve( df, score_cols = clf_list[[k]], lab_col = lab_col, fold_col = fold_col )	
		}
		else{
			attr = average_prc_curve( df, score_cols = clf_list[[k]], lab_col = lab_col, fold_col = fold_col )	
		}

		if( k == 1 )
		{						
			out_df = data.frame( x = attr[[1]]$x , attr[[1]]$y_avg  )
		}
		else
		{
			if( identical( out_df$x, attr[[1]]$x ) )
			{
				out_df = cbind ( out_df	, attr[[1]]$y_avg )		
			}
		}
	}
	colnames( out_df ) = c( "x", names(clf_list) )

	return(out_df)
}


## Sub sample ROC and PRC metric 
## Function :: sample_points
##	
##	Description : This function is sub-samples PRC or ROC metric
##

sample_points = function( x=NULL, y=NULL,  n = 1000)
{
	x_idx_first = 1
	x_idx_last = length(x);
	n = n - 2;
	sample_space = x_idx_first+1 : x_idx_last-1;
	random_idx = sample( sample_space, n, replace = F )
	random_sorted_idx = sort(random_idx);

	random_all_idx = c( x_idx_first, random_sorted_idx, x_idx_last )
	x_random_val = x[ random_all_idx ]
	y_random_val = y[ random_all_idx ];

	return( list( x_random_val = x_random_val, y_random_val = y_random_val ) );
}


## Call various performance performance function 
## Function :: call_performance_metric
##	
##	Description : This function is wrapper to call various performance metric function 
##

call_performance_metric = function( df = NULL, classifier = NULL, fold_column = "Fold", score_column = NULL , label_column = NULL, negative_source = NULL )
{
	p_at_recall_x = precision_at_recall_x( df = df, fold_column = fold_column, score_column = score_column, label_column = label_column )
	p_at_recall_x = data.frame( p_at_recall_x, Classifier = classifier, Negative_Source = negative_source )
	
	all_metric_list = prc_roc_auc(df = df, fold_column = fold_column, score_column = score_column, label_column = label_column )
	all_metric_df = data.frame( all_metric_list$prc_roc, Classifier = classifier , Negative_Source = negative_source )

	sampledPRC = metric_at_value_x( df = df, fold_column = fold_column, score_column = score_column, label_column = label_column, sampling_points = 20, metric = "PRC" )
	sampledPRC_df = data.frame( sampledPRC, Classifier = classifier , Negative_Source = negative_source )

	sampledROC = metric_at_value_x( df = df, fold_column = fold_column, score_column = score_column, label_column = label_column, sampling_points = 20, metric = "ROC" )
	sampledROC_df = data.frame( sampledROC, Classifier = classifier , Negative_Source = negative_source )

	metric_list = list( p_at_recall_x = p_at_recall_x, all_metric_df = all_metric_df, sampledPRC = sampledPRC_df, sampledROC = sampledROC_df );

	return(metric_list);
}


varying_neg_set = function( working_dir = NA , file = NA, label_column = NA, segment_column = NA, feature_start_col = NA, feature_end_col = NA, classifiers = c("rfc","svc"), neg_set_size = c( 1,2,5,10,20) , label = "Breast_Normal")
{
	setwd(working_dir);
	df = as.data.frame ( fread( file, sep = "\t", header = T, stringsAsFactors = F, check.names = F ) ) 

	pos_df	= df [ df[ , label_column ] == 2 , ] 
	neg_df	= df [ df[ , label_column ] == 0 , ] 

	neg_set_count = dim(pos_df)[1] * neg_set_size;

	for( i in 1:length( neg_set_count ) )
	{
		if( neg_set_count[i] > dim(neg_df)[1] )
		{
			print("Negative sub-sample size is larger than negative data set");
		}
		else
		{
			# Set the working directory #
			neg_set_dir = paste( working_dir , "/", neg_set_size[i] , "_times" , sep = "" );	
			dir.create( neg_set_dir )
			setwd(neg_set_dir)	

			neg_df_temp =  neg_df [  sample( dim(neg_df)[1], neg_set_count[i] , replace = F ) ,  ]
			merged_df = rbind( pos_df, neg_df_temp );
			merged_df = merged_df[ order( merged_df[,1], merged_df[,2] ) , ]
			merged_df[,label_column] = ifelse( merged_df[,label_column] == 2, 1, 0 )
			temp_file_name = paste( label, "_", neg_set_size[i] , "_times.txt" , sep = "" );
			write.table( merged_df , temp_file_name , sep = "\t", quote = F, row.names = F );

			for( m in 1:length( classifiers ) )
			{
				# Set the working directory #
				clf_dir = paste( neg_set_dir , "/", classifiers[m] , sep = "" );	
				dir.create( clf_dir )
				setwd(clf_dir)
				
				if( classifiers[m] == "rfc" )
				{
					file.symlink( paste( "../" , temp_file_name , sep = "" ) , temp_file_name  );

					rfc_command = paste( "bsub -q long -e " , neg_set_size[i] , ".e -o " , neg_set_size[i] , ".o -R \"select[type==X86_64 && mem > 20000] rusage[mem=20000]\" -M20000 '/software/team113/Python-3.4.3/bin/python3.4 /nfs/users/nfs_r/rm8/My_Code/Mutation_Annotation_Pipeline/Python/Python_Implementation/Mutation_Classification/segment_preserved_generic.py -i " ,temp_file_name , " -o ", getwd() , "/ -l ", label_column, " -s ", feature_start_col , " -e " , feature_end_col ," -f 5 -c ", segment_column ," -x rfc -w False -a 0.5 -b 0.5'", sep = "" );
					system(rfc_command)
				}
				else if( classifiers[m] == "svc" )
				{
					# ----- SVC ----- #
					#svc_command = paste( "bsub -q long -e " , neg_set_size[i] , ".e -o " , neg_set_size[i] , ".o -R \"select[type==X86_64 && mem > 20000] rusage[mem=20000]\" -M20000 '/software/team113/Python-3.4.3/bin/python3.4 /nfs/users/nfs_r/rm8/My_Code/Mutation_Annotation_Pipeline/Python/Python_Implementation/Mutation_Classification/segment_preserved.py -i " ,file , " -o ", getwd() , "/ -l Binom_Label -s 22 -e 314 -f 5 -c ", segment_column ," -x svc -w False -a 0.1 -b 0.9'", sep = "" );	
					#system(svc_command)	
				}
			
				# Reset Working Directory to Negative Set Size ##
				setwd(neg_set_dir)
			}
		
			# Reset Working Directory to main working directory ##
			setwd(working_dir)
		}
	}
}



## ============================================================================ ##
## 			 Function ::  Computes K-M Style Score for each feature
##	1. Not Rescaled 
## ============================================================================ ##


compute_running_score = function( df = NULL, feature = NULL, feature_threshold = 0, class_value = 1, class_label = "Binom_Label" )
{

	#print( feature );
	#print( class_value );

	df = df[ !is.na(rowSums( df ) ), ];

	if( length( unique( df[,class_label] ) ) == 1 )
	{	
		#print("All Feature value belong to one class")
		summary = list ( pos_class_index = 0 , neg_class = 0 , max_score = -1, feature = "NA", direction = "All_Same" );
		return( summary );
		
	}else if( length ( unique( df[,feature] ) ) == 2 && sum( unique( df[,feature] ) ) == 1 )
	{
		type = "Binary";

	}else if ( length ( unique( df[,feature] ) ) > 2 )
	{
		type = "Continuous";
	}
	else if( length ( unique( df[,feature] ) ) == 2 )
	{
		type = "Continuous";

	}
	else
	{
		#print("Can not determine feature type :: Feature will not be used")
		summary = list ( pos_class_index = 0 , neg_class = 0 , max_score = -1, feature = feature, direction = "Confused" )
		return( summary );
	}

	total_mutations = dim(df)[1]
	total_positive_mutations = dim( df[ df[,class_label] == class_value, ])[1]
	total_negative_mutations = dim( df[ df[,class_label] != class_value, ])[1]

	## ========================================================== ##
	##   Binary :: Subset Based on Feauture Threshold  = 0 or 1   ##

	if( type == "Binary" )
	{
		## -------------------- ##
		print("Running Binary Feature")
		
		df_o = df[order(df[,feature],decreasing=TRUE),]
		#print(df_o)
		forward_index = as.numeric( rownames( df_o[ df_o[,feature] != feature_threshold, ] ) )
		reverse_index = as.numeric( rownames( df_o[ df_o[,feature] == feature_threshold, ] ) )

		## Check if all the samples belong to one class ##
		if( length( unique( df[,class_label] ) ) == 1 )
		{	
			print("All Feature value belong to one class")
			summary = list ( pos_class_index = 0 , neg_class = 0 , max_score = -1, feature = feature, direction = "All_Same" );
		}
		else
		{
			print("There is impurity. Compute forwards and reverse score")


			#print(df_o[ df_o[,feature] != feature_threshold,  ])
			#print(df_o[ df_o[,feature] == feature_threshold,  ])

			forward_pos_sum = sum( df_o[ df_o[,feature] != feature_threshold, class_label ] )
			rev_pos_sum = sum( df_o[ df_o[,feature] == feature_threshold, class_label ] )

			#print(forward_pos_sum)
			#print(rev_pos_sum)

			forward_pos_mutation_score = sum( df_o[ df_o[,feature] != feature_threshold, class_label ] ) / total_positive_mutations	
			reverse_pos_mutation_score = sum( df_o[ df_o[,feature] == feature_threshold, class_label ] ) / total_negative_mutations
			
			print(total_positive_mutations)

			print(forward_pos_mutation_score)
			print(reverse_pos_mutation_score)

			if( forward_pos_mutation_score > reverse_pos_mutation_score )
			{	summary = list ( pos_class_index = forward_index , neg_class = reverse_index , max_score = forward_pos_mutation_score, feature = feature, direction = "forward" )
			}else
			{
				summary = list ( pos_class_index = reverse_index , neg_class = forward_index , max_score = reverse_pos_mutation_score, feature = feature, direction = "reverse")
			}
		}
		return( summary )
	}
	else
	{
		## ------------------------------ ##
		print("Running Continuous Feature")


		## Check if all the samples belong to one class ##
		if( length( unique( df[,class_label] ) ) == 1 )
		{	
			#print("All Feature value belong to one class")
			summary = list ( pos_class_index = 0 , neg_class = 0 , max_score = -1, feature = feature, direction = "All_Same" );
		}
		else
		{
		
			## ================= ##	
			## Forward Direction ##

			#print("Running Continuous Mode  :: Forward Direction ")

			forward_running_score = rep(0,total_mutations)		
			df_f = df[order(df[,feature],decreasing=TRUE),]
			print(df_f)

			for( i in 1:total_mutations )
			{
				forward_pos_mutation_score = 0 
				forward_neg_mutation_score = 0
				
				for( j in 1:i )
				{	
					if(  df_f[j,class_label] == class_value )
					{
						forward_pos_mutation_score = forward_pos_mutation_score + 1
					}
					else
					{
						forward_neg_mutation_score = forward_neg_mutation_score + 1
					}	
				}
				
				print( paste("Forward_Pos :", forward_pos_mutation_score, sep = "") )
				print( paste("Forward_Neg :", forward_neg_mutation_score, sep = "") )
				

				norm_forward_pos_mutation_score = forward_pos_mutation_score / total_positive_mutations
				norm_forward_neg_mutation_score = forward_neg_mutation_score / total_negative_mutations
				forward_running_score[i] = norm_forward_pos_mutation_score - norm_forward_neg_mutation_score

				print( paste("Forward Running Score :", forward_running_score[i], sep = "") )
			}

			## ================= ##	
			## Reverse Direction ##

			#print("Running Continuous Mode  :: Reverse Direction ")

			reverse_running_score = rep(0,total_mutations)		
			df_r = df[order(df[,feature],decreasing=FALSE),]
			print(df_r)

			for( i in 1:total_mutations )
			{
				reverse_pos_mutation_score = 0 
				reverse_neg_mutation_score = 0
				
				for( j in 1:i )
				{	
					if(  df_r[j,class_label] == class_value )
					{
						reverse_pos_mutation_score = reverse_pos_mutation_score + 1
					}
					else
					{
						reverse_neg_mutation_score = reverse_neg_mutation_score + 1
					}	
				}
				
				print( paste("Rev_Pos :", reverse_pos_mutation_score, sep = "") )
				print( paste("Rev_Neg :", reverse_neg_mutation_score, sep = "") )

				norm_reverse_pos_mutation_score = reverse_pos_mutation_score / total_positive_mutations
				norm_reverse_neg_mutation_score = reverse_neg_mutation_score / total_negative_mutations
				reverse_running_score[i] = norm_reverse_pos_mutation_score - norm_reverse_neg_mutation_score

				print( paste("Reverse Running Score :", reverse_running_score[i], sep = "") )
			}	

			## Compare Forward and Reverse Direction Score and Find Max ##
			if( max( forward_running_score ) >= max( reverse_running_score ) )
			{
				max_score = max( forward_running_score );
				forward_index = as.numeric( rownames( df_f[1:which.max( forward_running_score),] ) )
				forward_index = forward_index[ !is.na(forward_index) ]
				reverse_index = as.numeric( rownames( df_f[ which.max( forward_running_score ) + 1 : dim(df_f)[1] ,] ) )
				reverse_index = reverse_index[ !is.na(reverse_index) ]
				
				summary = list ( pos_class_index = forward_index , neg_class = reverse_index , max_score = max_score , feature = feature, direction = "forward" )
			}
			else
			{
				max_score = max( reverse_running_score );
				forward_index = as.numeric( rownames( df_r[1:which.max( reverse_running_score),] ) )
				forward_index = forward_index[ !is.na(forward_index) ]
				reverse_index = as.numeric( rownames( df_r[ which.max( reverse_running_score ) + 1 : dim(df_r)[1] ,] ) )
				reverse_index = reverse_index[ !is.na(reverse_index) ]

				summary = list ( pos_class_index = forward_index , neg_class = reverse_index , max_score = max_score, feature = feature, direction = "reverse" )
			}

			return ( summary );
		}
		
	} ## End of Continuous Else ##
}



## ============================================================================ ##
## 			 Function :: Performance Metric across different negative set
## ============================================================================ ##

performance_metric_var_neg = function( path = NULL, classifier = NULL, neg_set_size = NULL, score_column = NULL, label_column = NULL  )
{

	require(ROCR)

	clf_path = paste( path, "/", classifier, sep = "" )
	CLF_Files = list.files( path = clf_path , full.names= T, patter = ".txt" );

	print(CLF_Files)

	if( is.null(neg_set_size) )
	{
		neg_set_size = c("1K", "2K", "5K", "10K", "50K", "100K");	
	}
		
	tpr = c()
	fpr = c()
	precision = c()
	recall = c()
	label = c()
	recall_prec_sampled_list_df = list();
	tp_fp_sampled_list_df = list();

	for( i in 1:length( neg_set_size ) )
	{
		print(i);
		file = CLF_Files [ grep( paste( "_" , neg_set_size[i], ".", sep = "" ) , CLF_Files ) ]
		if( classifier == "NonConvex" )
		{
			df = read.delim( file, sep = "\t", header = F, stringsAsFactors = F, check.names = F )	
		}else{
			df = read.delim( file, sep = "\t", header = T, stringsAsFactors = F, check.names = F )	
		}

		if( classifier == "OCSVM" )
		{
			# Scale OCSVM Prediction Score and Swap Binom Label Back normal #
			df[,label_column] = df[,label_column] * -1 
    		
    		predicted_score = df[ , score_column ]
    		predicted_prob = rescale( predicted_score, 0 , 1 );	
    		predicted_prob_swapped = 1 - predicted_prob
    		df[ , score_column ] = predicted_prob_swapped
		}
		
		pred_df <- prediction(df[,score_column], df[,label_column])
		tp_fp_df <- performance(pred_df, "tpr", "fpr")
		pr_rc_df <- performance(pred_df, "prec", "rec")
		auc_df <- unlist ( slot( performance(pred_df, "auc") , "y.values" ) )
		prbe_df <- unlist ( slot( performance(pred_df, "prbe") , "y.values" ) )

		recall_df = unlist(slot(pr_rc_df, "x.values") )
		prec_df = unlist(slot(pr_rc_df, "y.values") )
		
		##   The first precision value is often  Inf | NaN   # 
		## So replace 1st precision value with Second value   ##
		
		prec_df[1] = prec_df[2]
		
		recall_prec_sampled_list_df[[i]] = vector_pair_subset( recall_df, prec_df , size = 200 );
		recall = append( recall, unlist ( recall_prec_sampled_list_df[[i]][ 1 ] ) )
		precision = append( precision, unlist ( recall_prec_sampled_list_df[[i]][ 2 ] ) )

		fpr_df = unlist(slot(tp_fp_df, "x.values") )
		tpr_df = unlist(slot(tp_fp_df, "y.values") )
		tp_fp_sampled_list_df[[i]] = vector_pair_subset( fpr_df, tpr_df , size = 200 );
		fpr = append( fpr , unlist( tp_fp_sampled_list_df[[i]][1] ) )
		tpr = append( tpr , unlist( tp_fp_sampled_list_df[[i]][2] ) )

		label = append( label, rep( neg_set_size[i] , 200 ) )
	}
	return( list( recall = recall, precision = precision, fpr = fpr, tpr = tpr , label = label ) );
}

## ====================================================================== ##
## 			 Function :: All Performance Metric 
## ====================================================================== ##


all_performance_metric_cancer = function( df = NULL, data_label = NULL, score_column = NULL, label_column = NULL, sampling_size = 20  )
{
	require(ROCR)

	tpr = c()
	fpr = c()
	precision = c()
	recall = c()
	label = c()
	
	pred_df <- prediction(df[,score_column], df[,label_column])
	
	tp_fp_df <- performance(pred_df, "tpr", "fpr")
	pr_rc_df <- performance(pred_df, "prec", "rec")
	
	auc_df <- unlist ( slot( performance(pred_df, "auc") , "y.values" ) )
	prbe_df <- unlist ( slot( performance(pred_df, "prbe") , "y.values" ) )

	mcc_x <- unlist ( slot( performance(pred_df, "mat") , "x.values" ) )
	mcc_y <- unlist ( slot( performance(pred_df, "mat") , "y.values" ) )
	mcc_x = mcc_x[ 2:(length(mcc_x) - 1) ]
	mcc_y = mcc_y[ 2:(length(mcc_y) - 1) ]	

	fscore_x <- unlist ( slot( performance(pred_df, "f") , "x.values" ) )	
	fscore_y <- unlist ( slot( performance(pred_df, "f") , "y.values" ) )	
	fscore_x = fscore_x[ 2:length(fscore_x) ]
	fscore_y = fscore_y[ 2:length(fscore_y) ]	
	
	cost_x <- unlist ( slot( performance(pred_df, "cost") , "x.values" ) )	
	cost_y <- unlist ( slot( performance(pred_df, "cost") , "y.values" ) )	
	cost_x = cost_x[ 2:length(cost_x) ]
	cost_y = cost_y[ 2:length(cost_y) ]	

	
	## Precision Recall 
		recall_df = unlist(slot(pr_rc_df, "x.values") )
		prec_df = unlist(slot(pr_rc_df, "y.values") )
		##   The first precision value is often  Inf | NaN   # 
		## So remove first recall and precision value    ##
		prec_df = prec_df[ 2:length(prec_df) ]
		recall_df = recall_df[ 2:length(recall_df) ]
		# --- Do not fix first and last value --- #
		# prec_df[1] = 1 ;  prec_df[ length(prec_df)] = 0;
		# recall_df[1] = 0 ; recall_df[ length(recall_df)] = 1;
		recall_prec_sampled_list_df = vector_pair_subset( recall_df, prec_df , size = sampling_size );
		recall = unlist ( recall_prec_sampled_list_df[[1]] )
		precision = unlist ( recall_prec_sampled_list_df[[2]] ) 

	## ROC 	
		fpr_df = unlist(slot(tp_fp_df, "x.values") )
		tpr_df = unlist(slot(tp_fp_df, "y.values") )
		tp_fp_sampled_list_df = vector_pair_subset( fpr_df, tpr_df , size = sampling_size );
		fpr = unlist( tp_fp_sampled_list_df[[1]] )
		tpr = unlist( tp_fp_sampled_list_df[[2]] ) 

	## MCC 	
		mcc_sampled_list_df = vector_pair_subset( mcc_x, mcc_y , size = sampling_size );
		mccX = unlist( mcc_sampled_list_df[[1]] )
		mccY = unlist( mcc_sampled_list_df[[2]] ) 

	## Fscore	
		fscore_sampled_list_df = vector_pair_subset( fscore_x, fscore_y , size = sampling_size );
		fscoreX = unlist( fscore_sampled_list_df[[1]] )
		fscoreY = unlist( fscore_sampled_list_df[[2]] ) 	

	## Cost	
		cost_sampled_list_df = vector_pair_subset( cost_x, cost_y , size = sampling_size );
		costX = unlist( cost_sampled_list_df[[1]] )
		costY = unlist( cost_sampled_list_df[[2]] ) 		

	label = rep( data_label , sampling_size )
	
	temp_df = data.frame( X = c( recall, fpr, mccX, fscoreX, costX ), Y = c( precision, tpr, mccY, fscoreY, costY ), Category = c( rep("PRC",sampling_size ), rep("ROC",sampling_size ), rep("MCC",sampling_size ), rep("Fscore",sampling_size ), rep("Cost",sampling_size ) ), Label = rep( data_label , sampling_size * 5 ) )

	#temp_df = data.frame( Recall = recall, Precision = precision, FPR = fpr, TPR = tpr, Label = label );
	return( temp_df );
}


## ====================================================================== ##
## 			 Function :: Basic Performance Metric 
## ====================================================================== ##


basic_performance_metric_cancer = function( df = NULL, data_label = NULL, score_column = NULL, label_column = NULL, sampling_size = 20  )
{
	require(ROCR)

	tpr = c()
	fpr = c()
	precision = c()
	recall = c()
	label = c()
	
	pred_df <- prediction(df[,score_column], df[,label_column])
	
	tp_fp_df <- performance(pred_df, "tpr", "fpr")
	pr_rc_df <- performance(pred_df, "prec", "rec")
	
	auc_df <- unlist ( slot( performance(pred_df, "auc") , "y.values" ) )
	prbe_df <- unlist ( slot( performance(pred_df, "prbe") , "y.values" ) )

	
	## Precision Recall 
		recall_df = unlist(slot(pr_rc_df, "x.values") )
		prec_df = unlist(slot(pr_rc_df, "y.values") )
		##   The first precision value is often  Inf | NaN   # 
		## So remove first recall and precision value    ##
		prec_df = prec_df[ 2:length(prec_df) ]
		recall_df = recall_df[ 2:length(recall_df) ]
		# --- Do not fix first and last value --- #
		# prec_df[1] = 1 ;  prec_df[ length(prec_df)] = 0;
		# recall_df[1] = 0 ; recall_df[ length(recall_df)] = 1;
		
		if( length(prec_df) < sampling_size )
		{
		 	prec_sampling = length(prec_df);	
		 }
		 else{
		 	prec_sampling = sampling_size
		}
		
		## With subsampling 	
		recall_prec_sampled_list_df = vector_pair_subset( recall_df, prec_df , size = prec_sampling );
		recall = unlist ( recall_prec_sampled_list_df[[1]] )
		precision = unlist ( recall_prec_sampled_list_df[[2]] ) 
		
		## Not doing subsampling 	
		# recall = recall_df;
		# precision = prec_df
		# prec_sampling = length(precision)

	## ROC 	
		fpr_df = unlist(slot(tp_fp_df, "x.values") )
		tpr_df = unlist(slot(tp_fp_df, "y.values") )

		
		if( length(fpr_df) < sampling_size )
		{
			fpr_sampling = length(fpr_df);
		}
		else{
			fpr_sampling = sampling_size
		}

		## With subsampling 	
		tp_fp_sampled_list_df = vector_pair_subset( fpr_df, tpr_df , size = fpr_sampling );
		fpr = unlist( tp_fp_sampled_list_df[[1]] )
		tpr = unlist( tp_fp_sampled_list_df[[2]] ) 

		## Without subsampling 	
		# fpr = fpr_df;
		# tpr = tpr_df;

		fpr_sampling = length(fpr)

		temp_df = data.frame( X = c( recall, fpr ), Y = c( precision, tpr ), Category = c( rep("PRC",prec_sampling ), rep("ROC",fpr_sampling ) ), Label = rep( data_label , fpr_sampling + prec_sampling ) )

	#temp_df = data.frame( Recall = recall, Precision = precision, FPR = fpr, TPR = tpr, Label = label );
	return( temp_df );
}


## ==================================================================================== ##
## 		Function :: Returns average performance metric of of N fold cross validation
##
##			Input :: 
##				1. Data frame with Actual Label and Prediction Score
##				2. Label for for the metric
##				3. Fold column
##				4. Score column
##				5. Actual label column
##				6. Sampling size
##
## ==================================================================================== ##

# df = scored_df
# folds = c(1,2,3,4,5)
# data_label = "100K_asymRFC30"
# score_column = "asym_30feat_RFC_score_inv"
# fold_col = "Fold"
# label_column = "bf_Binom_Label"
# sampling_size = 20


fold_mean_metric = function( df = NULL, folds = c(1,2,3,4,5), data_label = NULL, fold_col = NULL, score_column = NULL, label_column =  NULL, sampling_size = 20 )
{
	for( f in 1:length(folds) )
	{
		temp_df = df[ df[, fold_col] == (folds[f] - 1),  ]
		temp_met = basic_performance_metric_cancer( df = temp_df, data_label = data_label, score_column = score_column, label_column = label_column, sampling_size = sampling_size  )

		#a_mean = mean( mean(a_0) , mean(a_1) , mean(a_2) , mean(a_3) , mean(a_4) )
		#a_sd = sd( c( mean(a_0) , mean(a_1) , mean(a_2) , mean(a_3) , mean(a_4) ) )
		#a_se = a_sd / sqrt( 5 )

		if( f == 1 )
		{
			met = temp_met;
			PRE = temp_met[ temp_met$Category == "PRC", "Y" ]
			REC = temp_met[ temp_met$Category == "PRC", "X" ]

			TPR = temp_met[ temp_met$Category == "ROC", "Y" ]
			FPR = temp_met[ temp_met$Category == "ROC", "X" ]
		}
		else
		{
			met$X = met$X + temp_met$X
			met$Y = met$Y + temp_met$Y	

			PRE = cbind( PRE, temp_met[ temp_met$Category == "PRC", "Y" ] )
			REC = cbind( REC, temp_met[ temp_met$Category == "PRC", "X" ] )

			TPR = cbind( TPR, temp_met[ temp_met$Category == "ROC", "Y" ] )
			FPR = cbind( FPR, temp_met[ temp_met$Category == "ROC", "X" ] )

		}
	}

	met$X = met$X / length( folds )	
	met$Y = met$Y/ length( folds )
	met_df = data.frame( X = met$X, Y = met$Y, met[, c("Category","Label") ] );

	met_list = list( met_df = met_df, PRE = PRE, REC = REC, TPR = TPR, FPR = FPR );
	return(met_list);

}


## ==================================================================================== ##
## 		Function :: Subsamples performance metric based on Performance metric size
##
##			Input :: 
##				1. Data frame with X = Recall or FPR ; Y = Precision or TPR
##				2. Start window : vector specifying start of an interval
##				3. End window : vector specifying end of an interval
##
## ==================================================================================== ##

subsample_performance_met = function( submet = NULL, start_window = NULL, end_window = NULL  )
{
	require(dplyr)
	met_val = c()
	for( i in 1:length( start_window ) )
	{
		temp_val = mean(submet[ submet$X >= start_window[i] & submet$X < end_window[i], "Y" ]);	
		met_val = append( met_val, temp_val )
	}
	return(met_val);
}

## ==================================================================================== ##
## 		Function :: Compliles ROC and PRC from all folds
##
##			Input :: 
##				1. Data frame with Actual Label and Prediction Score
##				2. Label for for the metric
##				3. Fold column
##				4. Score column
##				5. Actual label column
##				6. Sampling size
##
## ==================================================================================== ##

# df = scored_df
# folds = c(1,2,3,4,5)
# data_label = "100K_asymRFC30"
# score_column = "asym_30feat_RFC_score_inv"
# fold_col = "Fold"
# label_column = "bf_Binom_Label"
# sampling_size = 20


fold_metric_compile = function( df = NULL, folds = c(1,2,3,4,5), data_label = NULL, fold_col = NULL, score_column = NULL, label_column =  NULL, sampling_size = 20 )
{
	require(dplyr)
	
	sampling_points = 1:100 * 0.05;
	sampling_window_start = sampling_points - 0.01
	sampling_window_end = sampling_points + 0.01

	## Get All PRC and ROC values ##
	for( f in 1:length(folds) )
	{
		temp_df = df[ df[, fold_col] == folds[f],  ]
		#print(dim(temp_df))
		# df = temp_df, data_label = data_label, score_column = score_column, label_column = label_column, sampling_size = sampling_size;

		## Add dummy score value compute PRC / ROC
		if( length(unique( temp_df[,score_column] ) ) == 1 ) 
		{
			temp_df[sample( dim(temp_df)[1] , 5, replace = F ),score_column] = 0.0001
		}
		#print(table( temp_df[,score_column] , temp_df[,label_column]))
		temp_met = basic_performance_metric_cancer( df = temp_df, data_label = data_label, score_column = score_column, label_column = label_column, sampling_size = sampling_size  )

		# PRC :: sub sampled precision
		temp_met_PRC = temp_met[ temp_met$Category == "PRC", ]
		sampled_PRE = subsample_performance_met( submet = temp_met_PRC, start_window = sampling_window_start, end_window = sampling_window_end  )	
		fold_recall_prec = data.frame( Recall = sampling_points, Precision = sampled_PRE, Category = "PRC", Label = data_label, Fold =  folds[f] ) 

		# ROC :: sub sampled True positive rate
		temp_met_ROC = temp_met[ temp_met$Category == "ROC", ]
		sampled_TPR = subsample_performance_met( submet = temp_met_ROC, start_window = sampling_window_start, end_window = sampling_window_end  )			
		fold_fpr_tpr = data.frame( FPR = sampling_points, TPR = sampled_TPR, Category = "ROC", Label = data_label, Fold =  folds[f] ) 

		if( f == 1 )
		{
			met = cbind( temp_met, Fold = rep((folds[f]), dim(temp_met)[1] ) );
			recall_prec = fold_recall_prec
			tpr_fpr = fold_fpr_tpr
		}
		else
		{
			met = rbind( met, cbind( temp_met, Fold = rep((folds[f]), dim(temp_met)[1] ) ) );
			recall_prec = rbind( recall_prec, fold_recall_prec )
			tpr_fpr = rbind( tpr_fpr, fold_fpr_tpr )
		}
	}

	run_metrics_list = list( met = met, recall_prec = recall_prec, tpr_fpr = tpr_fpr );
	return(run_metrics_list);
}


## ==================================================================================== ##
## 		Function :: Returns average, min and max performance metric of 
##					N fold cross validation
##
##			Input :: 
##				1. Data frame with Actual Label and Prediction Score
##				2. Label for for the metric
##				3. Fold column
##				4. Score column
##				5. Actual label column
##				6. Sampling size
##
## ==================================================================================== ##

# require(data.table)
# original_file = "/Users/rm8/Sanger/experiments/Mutation_Pattern/Classifier_Performance/New_ExonF_Data/Binomial_Distribution/Trained_Classifiers/Clf_for_Paper/COSMIC_2017/Final_RFC_SVC_Repeat_Run_bfBinom_Stratified_custom_norm/20K/Genome_Segment_Label_1K_segment_marked/Run1/rand_pos_neg.svc_rfc_Clf_Scored.txt"; 
# scored_file = "/Users/rm8/Sanger/experiments/Mutation_Pattern/Classifier_Performance/New_ExonF_Data/Binomial_Distribution/Trained_Classifiers/Clf_for_Paper/COSMIC_2017/Final_RFC_SVC_Repeat_Run_bfBinom_Stratified_custom_norm/20K/Genome_Segment_Label_1K_segment_marked/Run1/rand_pos_neg.svc_rfc_asymrfc_Scored.txt"; 
# original = as.data.frame( fread( original_file, sep = "\t", header = T, stringsAsFactors = F, check.names = F ) , check.names = F );
# scored = as.data.frame( fread( scored_file, sep = "\t", header = F, stringsAsFactors = F, check.names = F ) )
#
# if( identical( original$Bound_Motif, scored$V1 ) )
# {
# 	asymRFC_Score_Column = 309;
# 	df = data.frame( original, asymRFC_Pred_Score = scored[,asymRFC_Score_Column], check.names = F )
	
# 	# Inverse the asymRFC Score #
# 	df$inverse_asymRFC_Pred_Score = 1 - df$asymRFC_Pred_Score
# 	df$inverse_asymRFC_Pred_Label = ifelse( df$inverse_asymRFC_Pred_Score >= 0.5, 1, 0 )
# }

# df$svc_scaled_score = rescale( df$svc_Pred_Score, 0 , 1 );
# folds = c(1,2,3,4,5) - 1;
# fold_col = "Fold";
# score_column = "inverse_asymRFC_Pred_Score";
# label_column =  "bf_Binom_Label";
# sampling_size = 200;


fold_metric = function( df = NULL, folds = c(1,2,3,4,5), fold_col = NULL, score_column = NULL, label_column =  NULL, sampling_size = 20 )
{
	#print(folds)
	cutoff = rev ( 0:sampling_size * ( 1 / sampling_size ) );
	fold_mat = matrix( 0, length(cutoff), length(folds) );
	tpr_mat = fold_mat ; fpr_mat = fold_mat; 
	prec_mat = fold_mat ; recall_mat = fold_mat;

	for( f in 1:length(folds) )
	{
		temp_df = df[ df[, fold_col] == folds[f],  ]
		temp_met = basic_metric( temp_df[,score_column], temp_df[,label_column], cutoff = cutoff );

		tpr_mat[,f] = temp_met[,"tpr"];
		fpr_mat[,f] = temp_met[,"fpr"];
		prec_mat[,f] = temp_met[,"prec"];
		recall_mat[,f] = temp_met[,"recall"];
	}

	# find the mean, max and min for each cut off #
	min_max_met_df = data.frame( 
		
		ROC_mean_tpr = rowMeans(tpr_mat, na.rm = TRUE),
		ROC_max_tpr = rowMax(tpr_mat), 
		ROC_min_tpr = rowMins(tpr_mat ), 
		
		ROC_min_fpr = rowMins(fpr_mat), 
		ROC_max_fpr = rowMax(fpr_mat), 
		ROC_mean_fpr = rowMeans(fpr_mat, na.rm = TRUE), 
		
		PRC_min_prec = rowMins(prec_mat), 
		PRC_max_prec = rowMax(prec_mat), 
		PRC_mean_prec = rowMeans(prec_mat, na.rm = TRUE),	

		PRC_min_recall = rowMins(recall_mat), 
		PRC_max_recall = rowMax(recall_mat), 
		PRC_mean_recall = rowMeans(recall_mat, na.rm = TRUE),
		cutoff = cutoff );

	# find the mean, max and min for each cut off #
	print("In Mean SD section")
	mean_sd_met_df = data.frame(
		ROC_mean_tpr = rowMeans(tpr_mat, na.rm = TRUE),
		ROC_max_tpr = rowMeans(tpr_mat, na.rm = TRUE) + rowSD(tpr_mat ), 
		ROC_min_tpr = rowMeans(tpr_mat, na.rm = TRUE) - rowSD(tpr_mat ), 
		
		ROC_mean_fpr = rowMeans(fpr_mat, na.rm = TRUE), 
		ROC_max_fpr = rowMeans(fpr_mat, na.rm = TRUE) + rowSD(fpr_mat), 
		ROC_min_fpr = rowMeans(fpr_mat, na.rm = TRUE) - rowSD(fpr_mat), 
		
		PRC_mean_prec = rowMeans(prec_mat, na.rm = TRUE),
		PRC_max_prec = rowMeans(prec_mat, na.rm = TRUE) + rowSD(prec_mat), 
		PRC_min_prec = rowMeans(prec_mat, na.rm = TRUE) - rowSD(prec_mat), 

		PRC_mean_recall = rowMeans(recall_mat, na.rm = TRUE),
		PRC_max_recall = rowMeans(recall_mat, na.rm = TRUE) + rowSD(recall_mat), 
		PRC_min_recall = rowMeans(recall_mat, na.rm = TRUE) - rowSD(recall_mat), 

		cutoff = cutoff 
	);

	mean_half_sd_met_df = data.frame(
		ROC_mean_tpr = rowMeans(tpr_mat, na.rm = TRUE),
		ROC_max_tpr = rowMeans(tpr_mat, na.rm = TRUE) + (rowSD(tpr_mat ) / 2), 
		ROC_min_tpr = rowMeans(tpr_mat, na.rm = TRUE) - (rowSD(tpr_mat ) / 2) , 
		
		ROC_mean_fpr = rowMeans(fpr_mat, na.rm = TRUE), 
		ROC_max_fpr = rowMeans(fpr_mat, na.rm = TRUE) + ( rowSD(fpr_mat) / 2) , 
		ROC_min_fpr = rowMeans(fpr_mat, na.rm = TRUE) - ( rowSD(fpr_mat) / 2) , 
		
		PRC_mean_prec = rowMeans(prec_mat, na.rm = TRUE),
		PRC_max_prec = rowMeans(prec_mat, na.rm = TRUE) + ( rowSD(prec_mat) / 2 ) , 
		PRC_min_prec = rowMeans(prec_mat, na.rm = TRUE) - ( rowSD(prec_mat) / 2 ) , 

		PRC_mean_recall = rowMeans(recall_mat, na.rm = TRUE),
		PRC_max_recall = rowMeans(recall_mat, na.rm = TRUE) + ( rowSD(recall_mat) / 2 ), 
		PRC_min_recall = rowMeans(recall_mat, na.rm = TRUE) - ( rowSD(recall_mat) / 2 ) , 

		cutoff = cutoff 
	);
		
	#return(mean_sd_met_df);
	#return(min_max_met_df);
	return( mean_half_sd_met_df );
}

## -----------------------------------------  ##
##		Square Confusion Matrix
##	
##		Input  : predictions score and label 
##		Output : confusion matrix 
## -----------------------------------------  ##

squareTable <- function(x,y) 
{
    x <- factor(x)
    y <- factor(y)

    commonLevels <- sort(unique(c(levels(x), levels(y))))

    x <- factor(x, levels = commonLevels)
    y <- factor(y, levels = commonLevels)

    table(x,y)
}

## --------------------------------------------  ##
##		Basic Classifier Performance Metric
##
##		
##	
##		Input  : predictions, label, operating points
##		Output : Returns TPR, FPR , Precision and Recall
## --------------------------------------------  ##


basic_metric = function( predictions , labels, cutoff = NULL )
{ 

	pred_nd_label = data.frame( prediction = predictions, labels = labels)

	tpr = c() ; fpr = c();
	prec = c() ; recall = c();

	for( i in 1:length( cutoff ) )
	{
		pred = ifelse( pred_nd_label$prediction >= cutoff[i], 1 , 0 )
		#print(pred)
		#print(labels)
		conf_mat = squareTable( labels, pred );
		#print(conf_mat)

		tp = conf_mat[2,2] ; fn = conf_mat[2,1] ;
		fp = conf_mat[1,2] ; tn = conf_mat[1,1] ;

		tpr = append( tpr, tp / (tp + fn) );
		fpr = append( fpr, fp / (fp + tn) );
		prec = append( prec, tp / (tp + fp ) );
		recall = append( recall, tp / (tp + fn) ); # Same as TPR
	}

	df = data.frame( tpr = tpr, fpr = fpr, prec = prec, recall = recall )
	return(df);
}


## ====================================================================== ##
## 			 Function :: Basic Performance Metric 
## ====================================================================== ##


precision_at_recall_x = function( df = NULL, fold_column = NULL, score_column = NULL, label_column = NULL )
{
	require(ROCR)

	#df = scored_df; fold_column = "Fold"; score_column = "rfc_gini_30feat_Pred_Score"; label_column = "bf_Binom_Label";	
	uniq_folds = unique( df[,fold_column] )
	rec_points = 6

	for( f in 1:length(uniq_folds) )
	{	
		df_fold = df[ df[,fold_column] == uniq_folds[f], ]

		#print( df_fold[1:10, c( label_column, score_column ) ] )

		pred_df <- prediction(df_fold[,score_column], df_fold[,label_column])
		pr_rc_df <- performance(pred_df, "prec", "rec")

		recall = unlist(slot(pr_rc_df, "x.values") )
		prec = unlist(slot(pr_rc_df, "y.values") )

		prec_005 = mean( prec[ recall > 0.03 & recall < 0.08 ] )
		prec_010 = mean( prec[ recall > 0.09 & recall < 0.11 ] )
		prec_015 = mean( prec[ recall > 0.12 & recall < 0.17 ] )
		prec_025 = mean( prec[ recall > 0.23 & recall < 0.27 ] )
		prec_050 = mean( prec[ recall > 0.45 & recall < 0.55 ] )
		prec_075 = mean( prec[ recall > 0.70 & recall < 0.80 ] )

		temp_prec_df = data.frame( fold = rep( uniq_folds[f], rec_points ), Recall = c( 0.05, 0.10, 0.15, 0.25, 0.50, 0.75 ), Precision = c( prec_005, prec_010, prec_015, prec_025, prec_050, prec_075 ) )

		if( f == 1)
		{
			precision_at_x_df = temp_prec_df
		}
		else{
			precision_at_x_df = rbind( precision_at_x_df, temp_prec_df )
		}
	}

	return(precision_at_x_df);
}


## ====================================================================== ##
## 			 Function :: Metric at value X
## ====================================================================== ##


metric_at_value_x = function( df = NULL, fold_column = NULL, score_column = NULL, label_column = NULL, sampling_points = sampling_points, metric = "PRC" )
{

	# df = scored_df; fold_column = "Fold"; score_column = "rfc_gini_30feat_Pred_Score"; label_column = "bf_Binom_Label"; sampling_points = 100; metric = "ROC"
	
	require(ROCR)
	uniq_folds = unique( df[,fold_column] )

	mid_point = ( 1:sampling_points ) * (1/sampling_points)
	lower_end = mid_point - 0.005
	upper_end = mid_point + 0.005
	
	for( f in 1:length(uniq_folds) )
	{	
		df_fold = df[ df[,fold_column] == uniq_folds[f], ]
		#print( df_fold[1:10, c( label_column, score_column ) ] )

		pred <- prediction(df_fold[,score_column], df_fold[,label_column])
		if( metric == "PRC" )
		{
			print("PRC")
			perform <- performance(pred, "prec", "rec")
			X = unlist( slot( perform, "x.values") )
			Y = unlist( slot( perform, "y.values") )	
		}else if( metric == "ROC" )
		{
			print("ROC")
			perform <- performance(pred, "tpr", "fpr")
			X = unlist( slot( perform, "x.values") )
			Y = unlist( slot( perform, "y.values") )	
		}
		
		## Define empty performance metric
		#print(X)
		#print(Y)
		#print(lower_end)
		#print(upper_end)
		mean_metric = c()
		for( j in 1:length( mid_point ) )
		{
			mean_metric = append( mean_metric, mean( Y[ X >= lower_end[j] & X <= upper_end[j] ] ) )
		}
				
		fold_metric = data.frame( fold = rep( uniq_folds[f], sampling_points ), X = mid_point, Y = mean_metric )

		if( f == 1)
		{
			metric_at_x_df = fold_metric
		}
		else{
			metric_at_x_df = rbind( metric_at_x_df, fold_metric )
		}
	}

	return(metric_at_x_df);
}



## ====================================================================== ##
## 			 Function :: Precision Recall, AUCPRC, AUCROC 
##
##	Intput :: 
##
##	1. dataframe with label and prediction score column
##	2. fold column
##	3. score column
##	4. label clumn
##
##	Output :: 
##		A list containing :
##	1. Precision Recall [ data frame ]
##	2. AUCPRC
##	3. AUCROC
##		
## ====================================================================== ##

#df = scored_df; fold_column = "Fold"; score_column = "rfc_gini_30feat_Pred_Score"; label_column = "bf_Binom_Label";
#df = scored_df; fold_column = "Fold"; score_column = "asym_30feat_RFC_score_inv"; label_column = "bf_Binom_Label";

prc_roc_auc = function( df = NULL, fold_column = NULL, score_column = NULL, label_column = NULL )
{
	require(ROCR)
	require(PRROC)

	uniq_folds = sort(unique( df[,fold_column] ))
	auc_prcurve = c();
	auc_roccurve = c();

	for( f in 1:length(uniq_folds) )
	{
		df_fold = df[ df[,fold_column] == uniq_folds[f], ]

		## Adding tiny varity in score to compute PRC / ROC
		if( length(unique( df_fold[,score_column] ) ) == 1 ) 
		{
			df_fold[sample( dim(df_fold)[1] , 5, replace = F ),score_column] = 0.0001
		}

		pos_score = df_fold[ df_fold[,label_column] == 1, score_column ]  
		neg_score = df_fold[ df_fold[,label_column] == 0, score_column ]  
		pr_curve_stats = pr.curve(scores.class0 = pos_score, scores.class1 = neg_score, curve=T)
		auc_prcurve = append( auc_prcurve, pr_curve_stats$auc.integral )

		roc_curve_stats = roc.curve(scores.class0 = pos_score, scores.class1 = neg_score, curve=T)
		auc_roccurve = append( auc_roccurve, roc_curve_stats$auc )


		pred_df <- prediction(df_fold[,score_column], df_fold[,label_column])
		
		# precision recall 
		pr_rc_df <- performance(pred_df, "prec", "rec")
		recall = unlist(slot(pr_rc_df, "x.values") )
		prec = unlist(slot(pr_rc_df, "y.values") )

		# precision recall 
		tpr_fpr_df <- performance(pred_df, "tpr", "fpr")
		fpr = unlist(slot(tpr_fpr_df, "x.values") )
		tpr = unlist(slot(tpr_fpr_df, "y.values") )

		# f-measure
		f_measure = performance(pred_df, "prbe" )


		cost_x <- unlist ( slot( performance(pred_df, "cost") , "x.values" ) )	
		cost_y <- unlist ( slot( performance(pred_df, "cost") , "y.values" ) )	
		cost_x = cost_x[ 2:length(cost_x) ]
		cost_y = cost_y[ 2:length(cost_y) ]

		temp_prc_df = data.frame( fold = rep( uniq_folds[f], length(recall) ), X = recall, Y = prec, Label = "PRC" )
		temp_roc_df = data.frame( fold = rep( uniq_folds[f], length(fpr) ), X = fpr, Y = tpr, Label = "ROC" )

		temp_roc_prc = rbind( temp_prc_df, temp_roc_df )

		if( f == 1)
		{
			prc_roc = temp_roc_prc
		}
		else{
			prc_roc = rbind( prc_roc, temp_roc_prc )
		}
	}
	pr_met = list( prc_roc = prc_roc, auc_prcurve = auc_prcurve, auc_roccurve = auc_roccurve );
	return(pr_met);
}

