#  Example of use below function
remove.NA.rows<-function(x) {
	NAs.per.row<-apply(x,1,sum.nas<-function(y) {sum(is.na(y))})
	return(x[NAs.per.row==0,])
}

#  summary(airquality)  #  note NA's
#  airquality.no.na<-remove.NA.rows(airquality)
#  summary(airquality.no.na)  #  note no NA's

