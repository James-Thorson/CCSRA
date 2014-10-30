

### CleanADMB
CleanAdmbFn=function(SaveFile=NULL, PrefixVec, KeepVec=c("pin","par","dat","rep"),quiet=FALSE){
  if(is.null(SaveFile)){SaveFile = getwd()}
  Index = 0

  for(i in 1:length(PrefixVec)){
    List=list.files(SaveFile)
    Grep = grep(PrefixVec[i],List)
    if(length(Grep)>0){
      Remove = List[Grep]
      Keep = vector(length=0)
      for(j in 1:length(KeepVec)){
        Keep = union(Keep,List[grep(KeepVec[j],List)])
      }
      Remove = setdiff(Remove,Keep)
      if(length(Remove)>0){
        file.remove(paste(SaveFile,Remove,sep=""))
      }
      Index = Index+length(Remove)
    }
  }
  Remove = c("eigv.rpt","fmin.log","variance","cmpdiff.tmp","gradfil1.tmp","admodel.cov","admodel.hes","admodel.dep")
    Remove = Remove[which(Remove %in% list.files(SaveFile))]
    if(length(Remove)>0) file.remove(paste(SaveFile,Remove,sep=""))
  Remove = c("f1b2list1","f1b2list12","f1b2list13","nf1b2list1","nf1b2list12","nf1b2list13")
    Remove = Remove[which(Remove %in% list.files(SaveFile))]
    if(length(Remove)>0) file.remove(paste(SaveFile,Remove,sep=""))
  if(quiet==FALSE) print(paste("Deleted ",Index," files from ",SaveFile,sep=""))
}

