df = read.delim("dummy_regression_data.csv", sep = ",")
ggplot(radial,aes(y=NTAV,x=age,color=factor(DM)))+geom_point()+stat_smooth(method="lm",se=FALSE)

