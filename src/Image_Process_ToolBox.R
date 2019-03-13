## ====================================================================== ##
##  Function to convert genotype data to binary values
##

genotype_binary_conversion = function( genotype )
{
	return( ifelse( genotype == "./." | genotype == "0/0", 0, 1 ) )
}


## ====================================================================== ##
##  Function to 
##  1. Rescale  
##  2. Add label 
##  3. Plot smoothed data
##


rescale_and_plot = function( labeled_data )
{
	print("Image Smoothing")	
	require(fields)
	X_scaled = rescale(labeled_data[,1], bottom = 0, top = 100)
	Y_scaled = rescale(labeled_data[,2], bottom = 0, top = 100)
	if( dim(labeled_data)[2] >= 4 )
	{
		feature_flag = TRUE	
		labled_scaled = data.frame(X=round(X_scaled), Y=round(Y_scaled), label = labeled_data[,3], feature =  labeled_data[,4] )
	}
	else{
		feature_flag = FALSE
		labled_scaled = data.frame(X=round(X_scaled), Y=round(Y_scaled), label = labeled_data[,3])	
	}
		
	pos_image_mat = matrix(0,100,100)
	neg_image_mat = matrix(0,100,100)

	for(i in 1:100)
	{
		for( j in 1:100 )
		{
			temp = labled_scaled[ labled_scaled[,1] == i & labled_scaled[,2] == j , ];
			## If any values overlap with this pixel ##
			if( dim(temp)[1] > 0 )
			{

				## Create A Pos Matrix ##
				pos_count = temp[ temp[,"label"] > 0,  ]
				if( dim( pos_count )[1] > 0 )
				{
					if (feature_flag)
					{
						pos_feature_count = pos_count [ pos_count[,"feature"] > 0 , ]
						if( dim( pos_feature_count )[1] > 0 )
						{
							pos_image_mat[i,j] = dim( pos_feature_count )[1]	
						}
						else
						{
							pos_image_mat[i,j] = 0
						}	
					}
					else
					{
						pos_image_mat[i,j] = dim( pos_count )[1]	
					}
				}	
				else
				{
					pos_image_mat[i,j] = 0	
				}
				
				## Create A Neg Matrix ##
				neg_count = temp[ temp[,"label"] == 0,  ]
				if( dim( neg_count )[1] > 0 )
				{
					if (feature_flag)
					{
						neg_feature_count = neg_count [ neg_count[,"feature"] > 0 , ]
						if( dim( neg_feature_count )[1] > 0 )
						{
							neg_image_mat[i,j] = dim( neg_feature_count )[1]	
						}
						else
						{
							neg_image_mat[i,j] = 0
						}		
					}
					else
					{
						neg_image_mat[i,j] = dim(neg_count)[1]
					}
				}	
				else
				{
					neg_image_mat[i,j] = 0	
				}
			}
			else
			{
				neg_image_mat[i,j] = 0	
				pos_image_mat[i,j] = 0	
			}
		}
	}	
	neg_image_mat_bk = neg_image_mat
	pos_image_mat_bk = pos_image_mat

	pos_norm = pos_image_mat / sum(pos_image_mat)
	neg_norm = neg_image_mat / sum(neg_image_mat)

	pos_norm = pos_norm + 1
	neg_norm = neg_norm + 1

	wght <- setup.image.smooth( nrow=100, ncol=100, dx=0.50, dy=0.50, theta=1, xwidth=0, ywidth=0 )
	pos_norm_smoothed <- image.smooth( pos_norm, dx=0.5, dy=0.5, wght)			
	neg_norm_smoothed <- image.smooth( neg_norm, dx=0.5, dy=0.5, wght)			

	image_ratio = pos_norm_smoothed$z / neg_norm_smoothed$z

	image_list = list( pos_image = pos_image_mat_bk, neg_image = neg_image_mat_bk, pos_image_smoothed = pos_norm_smoothed$z, neg_image_smoothed = neg_norm_smoothed$z, ratio_smoothed = image_ratio )
	return(image_list)
}

# look  <- image.smooth( test.image , xwidth=100, ywidth=100, dx=dx, dy=dx, theta=1.5 )
# look2 <- image.smooth( test.image , xwidth=1.5, ywidth=1.5, dx=dx, dy=dy, theta=1.5 )

# # compare these two
# set.panel( 1,2)
# image.plot( look)
# title("free boundaries")
# image.plot( look2) # look for periodic continuity at edges!
# title("periodic boundary in horizontal");
# set.panel(1,1);




rescale_and_plot_2d = function( labeled_data, sigma = 2 )
{
	#labeled_data = tSNE_clustered_annotated[,c(V1,V2,label)]
	print("Image Smoothing")	
	require(fields)
	require(EBImage)
	X_scaled = rescale(labeled_data[,1], bottom = 0, top = 100)
	Y_scaled = rescale(labeled_data[,2], bottom = 0, top = 100)
	if( dim(labeled_data)[2] >= 4 )
	{
		feature_flag = TRUE	
		labled_scaled = data.frame(X=round(X_scaled), Y=round(Y_scaled), label = labeled_data[,3], feature =  labeled_data[,4] )
	}else{
		feature_flag = FALSE
		labled_scaled = data.frame(X=round(X_scaled), Y=round(Y_scaled), label = labeled_data[,3])	
	}
	
	pos_image_mat = matrix(0,100,100)
	neg_image_mat = matrix(0,100,100)

	for(i in 1:100)
	{
		for( j in 1:100 )
		{
			temp = labled_scaled[ labled_scaled[,1] == i & labled_scaled[,2] == j , ];
			## If any values overlap with this pixel ##
			if( dim(temp)[1] > 0 )
			{

				## Create A Pos Matrix ##
				pos_count = temp[ temp[,"label"] > 0,  ]
				if( dim( pos_count )[1] > 0 )
				{
					if (feature_flag)
					{
						pos_feature_count = pos_count [ pos_count[,"feature"] > 0 , ]
						if( dim( pos_feature_count )[1] > 0 )
						{
							pos_image_mat[i,j] = dim( pos_feature_count )[1]	
						}
						else
						{
							pos_image_mat[i,j] = 0
						}	
					}
					else
					{
						pos_image_mat[i,j] = dim( pos_count )[1]	
					}
				}	
				else
				{
					pos_image_mat[i,j] = 0	
				}
				
				## Create A Neg Matrix ##
				neg_count = temp[ temp[,"label"] == 0,  ]
				if( dim( neg_count )[1] > 0 )
				{
					if (feature_flag)
					{
						neg_feature_count = neg_count [ neg_count[,"feature"] > 0 , ]
						if( dim( neg_feature_count )[1] > 0 )
						{
							neg_image_mat[i,j] = dim( neg_feature_count )[1]	
						}
						else
						{
							neg_image_mat[i,j] = 0
						}		
					}
					else
					{
						neg_image_mat[i,j] = dim(neg_count)[1]
					}
				}	
				else
				{
					neg_image_mat[i,j] = 0	
				}
			}
			else
			{
				neg_image_mat[i,j] = 0	
				pos_image_mat[i,j] = 0	
			}
		}
	}	
	
	neg_image_mat_bk = neg_image_mat
	pos_image_mat_bk = pos_image_mat

	pos_norm = pos_image_mat / sum(pos_image_mat)
	neg_norm = neg_image_mat / sum(neg_image_mat)

	pos_norm = pos_norm + 1
	neg_norm = neg_norm + 1

	pos_norm_gaus = gblur(pos_norm, sigma = sigma )
	neg_norm_gaus = gblur(neg_norm, sigma = sigma )
	
	image_ratio = pos_norm_gaus / neg_norm_gaus
	
	image_list = list( pos_image = pos_image_mat_bk, neg_image = neg_image_mat_bk, pos_image_smoothed = pos_norm_gaus, neg_image_smoothed = neg_norm_gaus, ratio_smoothed = image_ratio )
	return(image_list)
}

 