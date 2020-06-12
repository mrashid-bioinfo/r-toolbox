 ################################################
 #                                 				#
 #######    R Stat Process Toolbox    ###########
 #                                 				#
 ###############################################

CV <- function(x)
{
	  avg = mean(x)
	  x_std = sd(x)	
      return(x_std/avg)
}

z_normalize <- function(x)
{
	  mean_subtracted = x - mean(x)
	  x_std = sd(x)	
      return(mean_subtracted/x_std)
}

## Function to normalize a dataframe whose first column is feature_id

z_normalize_feat_df <- function(df=NULL)
{
	data = df[ , 2:dim(df)[2] ]
	data_z = apply( data, 2, function( x ){
		mean_subtracted = x - mean(x)
	  	x_std = sd(x)	
      return(mean_subtracted/x_std)	
	} )

	df_z = data.frame( feature_name = df[,1], data_z )
}


correlated_data_example = function()
{
	df = data.frame( a = sort( sample(1:50, 20,replace=T ) ), b= sort( sample(1:50,20,replace=T ) ), c = sort( sample(1:50,20,replace=T ), decreasing = T ), d = sample(1:50,20,replace=T ), e = sample(1:50,20,replace=T ), res = sort( sample(1:100, 20,replace=T ) ) )

	return(df)
}

 normalize_within_sample <- function( Matrix = NULL , method = c("standardise") )
 {	 
	 if( !is.null( Matrix ) )
	 {
		if( method == "standardise" ){
			
			#Matrix = 
			
		}	 
		 
	 }else{
		 print(" Matrix can not be NULL ");
		 return(0);
	 }
 } 



compute_delta = function( val = c(0,1,2,3) )
{
	
	delta = c()
	for( i in 1:length( val ) )
	{
		temp_val = val[-i]
		delta = append( delta, temp_val - val[i] )
	}
	return( delta [ !duplicated( abs( delta) ) ] )
}


######## Function getCIs
GetCIs <- function(Lab,Response,alpha=0.10)
{
	alpha1 = alpha/2 #alpha1 is the error rate for a one-sided confidence interval
	L = length(unique(Lab)) # L = No. of labs
	N = length(Lab) # N = Total number of tests across the labs
	lab.K = tapply(Response,Lab,length) #K in each lab
	#print(lab.K)
	lab.mean = tapply(Response,Lab,mean) #mean for each lab
	ystar = mean(lab.mean) # mean of lab means estimates overall mean
	KH = 1/mean(1/lab.K) # harmonic mean of K per lab
	dev = lab.mean-ystar #deviations of lab means from ystar
	#cat(paste0("Dev : ", dev , "\n") )
	MSU = KH*sum(dev^2)/(L-1) #BQI equation 11: MSamong using Unweighted SSamong
	#
	#ystar approx CI BQI (12): Overall Mean
	L12 = ystar - qt(1-alpha1,L-1)*sqrt(MSU/(L*KH))
	U12 = ystar + qt(1-alpha1,L-1)*sqrt(MSU/(L*KH))
	#
	#sr exact CI BQI (4): Repeatability SD
	#Use ANOVA (MOM) MSE, not REML MSE
	lab.var = tapply(Response,Lab,var)
	#cat(paste0("lab.var : ", lab.var , "\n") )
	lab.sd = sqrt(lab.var)
	#cat(paste0("lab.sd : ", lab.sd , "\n") )
	Kminus1 = lab.K - 1
	MSE = sum(lab.var*Kminus1)/sum(Kminus1)
	sr = sqrt(MSE)
	L4 = sqrt(MSE/F1inf(1-alpha1,N-L))
	U4 = sqrt(MSE/F1inf(alpha1,N-L))
	#
	#sR approx CI BQI (14): Reproducibility SD
	sR = sqrt((MSU/KH)+((KH-1)*MSE/KH)) # reproducibility SD using MSU
	G1 = 1-Finf2(alpha1,L-1)
	G2 = 1-Finf2(alpha1,N-L)
	H1 = Finf2(1-alpha1,L-1)-1
	H2 = Finf2(1-alpha1,N-L)-1
	L14 = sqrt(sR^2 - (1/KH)*sqrt(((G1^2)*MSU^2)+((G2^2)*((KH-1)^2)*MSE^2)))
	U14 = sqrt(sR^2 + (1/KH)*sqrt(((H1^2)*MSU^2)+((H2^2)*((KH-1)^2)*MSE^2)))
	#
	#r Wald CI BQI (15): intraclass correlation coefficient
	varA = max(0,(MSU - MSE)/KH)
	r = varA/(sR^2)
	L15 = (MSU/(KH*MSE*qf(1-alpha1,L-1,N-L)))-(1/min(lab.K))
	L15 = L15/(1+L15)
	L15 = max(0,L15)
	U15 = (MSU/(KH*MSE*qf(alpha1,L-1,N-L)))-(1/max(lab.K))
	U15 = U15/(1+U15)
	### end of calculations
	#
	### Provide output:
	d.collab = cbind(Lab,Response)
	colnames(d.collab)=c("Lab","Response")
	rownames(d.collab) = 1:length(Lab)
	rc = dim(d.collab)
	#
	### Print data checks and analysis results
	#cat("Number of labs and the total number of tests:","\n")
	#cat(c(L,N), "\n")
	#cat("\n")
	#
	##cat("First 3 and last 3 lines of data:", "\n")
	st = as.integer(rc[1] - 2)
	sf = as.integer(rc[1])
	#print(d.collab[c(1:3,st:sf),])
	##cat("\n")
	#
	#cat("Calculations performed with significance level alpha = ","\n")
	#cat(alpha,"\n")
	#cat("\n")
	#
	#cat("Number of tests in each lab:","\n")
	#print(lab.K)
	#cat("\n")
	#
	#cat("Harmonic mean of the number of tests at each lab:","\n")
	#cat(KH,"\n")
	#cat("\n")
	#
	#cat("Lab means:","\n")
	#print(lab.mean)
	#cat("\n")
	#
	#cat("Overall mean (=mean of lab means):","\n")
	#cat(ystar,"\n")
	#cat("\n")
	#
	#cat("Lab repeatability SDs","\n")
	#print(lab.sd)
	#cat("\n")
	#cat("Overall repeatability SD (pooled over labs)","\n")
	#cat(sr,"\n")
	#cat("\n")
	#
	#cat("Intermediate ANOVA calculations:","\n")
	#cat("Mean Square Among Labs, Mean Square within labs, Var among labs, interlab correlation","\n")
	#cat(c(MSU, MSE, varA, r),"\n")
	#cat("\n")
	#
	#cat("Intermediate calculations for approximate modified large sample method in BQI:","\n")
	#cat("G1, G2, H1, H2","\n")
	#cat(c(G1,G2,H1,H2),"\n")
	#cat("\n")
	#
	#cat("CONFIDENCE INTERVALS","\n")
	#cat("\n")
	#
	#cat("Overall mean response: Estimate, Lo, Up:","\n")
	#cat(c(ystar, L12, U12),"\n")
	#cat("\n")
	#
	cat("Repeatability SD: Estimate, Lo, Up:","\n")
	cat(c(sr,L4,U4),"\n")
	cat("\n")
	#
	#cat("Reproducibility SD: Estimate, Lo, Up:","\n")
	#cat(c(sR,L14,U14),"\n")
	#cat("\n")

	#
	#cat("Intralab correlation coefficient: Estimate, Lo, Up:","\n")
	#cat(c(r,L15,U15),"\n")
	#cat("\n")

	stats = list();
	stats[[ "Interlab_corr" ]] = c(r,L15,U15)
	stats[[ "Reproducibility" ]] = c(sR,L14,U14)
	stats[[ "Repeatability" ]] = c(sr,L4,U4)

	return( stats )

}
####### END OF GetCIs

######## NECESSARY FUNCTIONS
# These two quantile look-up functions are used in the formulas:
## F1inf ##
F1inf = function(p,df1)
{
	qc2 = qchisq(p,df1)/df1
	qc2
} # end of F1inf

## Finf2 ##
Finf2 = function(p,df2)
{
	qc1 = df2/qchisq(1-p,df2)
	qc1
} # end of Finf2
####################



performanace_analytics = function()
{
	require("PerformanceAnalytics")
	my_data <- mtcars[, c(1,3,4,5,6,7)]
	chart.Correlation(my_data, histogram=TRUE, pch=19)	
}
    


upper.panel<-function(x, y){
  points(x,y, pch=19, col=c("red", "green3", "blue")[iris$Species])
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  text(0.5, 0.9, txt)
}
