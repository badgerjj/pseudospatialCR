########################################################################
#FUNCTION TO BUILD CAPTURE HISTORY
captureHistories=function(dat=NULL, id.col=NULL, n.occ.col=NULL)
{
  #n.occ.col is the column specifying the unit to be used as an occasion (e.g. year)
  #determining the number of occasions
  occ=sort(unique(dat[, n.occ.col]),decreasing=FALSE) #creating a vector of each occasion
  n.occ=length(occ) #determining the number of occasions in the file

  ID=unique(dat[,id.col]) #vector with unique animal ID
  out=NULL #temporary output object

  for (i in 1:n.occ)
  {
    temp.dat=subset(dat, dat[,n.occ.col]==occ[i])
    temp.hist=ID%in%temp.dat[,id.col]
    temp.hist[FALSE]=0
    out=cbind(out, temp.hist)
  }

  ch=(apply(out,1,paste,collapse=""))
  ch=as.character(ch)
  output=data.frame(ID, ch, out, stringsAsFactors=FALSE)
  names(output)[3:(n.occ+2)]=as.character(occ)
  output=output[with(output, order(ID)),]

  #attaching Group covariate (Cluster in the case of the Baird's Pc analysis)
  cl=data.frame(dat$IndID, dat$Cluster)
  names(cl)=c("ID", "Cluster")
  output=merge(output, cl, by.x="ID", by.y="ID")
  final.output=unique(output)

  return(final.output)

}