
a = c(0.001,5,10)
b = c(0.001,200,500,1000)

for( i in 1:1000 )
{
	if ( sample( c("a","b"), 1, replace = F ) == "a" )
	{
		temp = sample(a,10,replace=T)
	}else{
		temp = sample(b,10,replace=T)
	}

	if( i == 1 )
	{
		ab = temp
	}
	else{
		ab = rbind( ab, temp )
	}
}

colnames(ab) = c("a1","a2","a3","a4","a5","b1","b2","b3","b4","b5")

cor_ab = cor(ab)
cor_log2_ab = cor(log2(ab))

par( mfcol=c(1,2) )
pheatmap(cor_ab, main="Raw" )
pheatmap(cor_log2_ab, main="Log2" )

dev.off()
