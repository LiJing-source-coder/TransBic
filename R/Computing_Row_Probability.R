#'
#'Perform the step of computing row probability p for all the numbers of 1 under the identified columns for a corresponding submatrix in the ranking matrix
#'
#'@param conditionList the condition partition
#'
#'@return colAndpr: the final identified columns of the corresponding submatrix (its number is denoted by Num(Col)) and all probability p for all the numbers of 1. To meet our destination, we merge them and the result is shown as follows:
#'the number of 1 (denoted by eNum): 0 and Num(Col), 1 and Num(Col)-1, 2 and Num(Col)-2, ..., Num(Col)/2(or(Num(Col)-1)/2 and (Num(Col)+1)/2)
#'the probability p:              P(0)+P(Num(Col)),P(1)+P(Num(Col)-1),P(2)+P(Num(Col)-2),...,P(Num(Col)/2)(or P((Num(Col)-1)/2)+P((Num(Col)+1)/2))
#'only output the probability p,not the number eNum.
#'
#'@export
#'
#'@examples
#'data(sim)
#'aa<-sample(col(sim),3)
#'bb<-sample(col(sim),2)
#'conList<-list(aa,bb)
#'prow_compute(conList)

"prow_compute"<-function(conditionList)
{
  fun_sig<-function(prenum,k)
  {
    fraList<-list(sapply((0:prenum),function(n) list(n)))
    fraList<-c(fraList,list(sapply((0:prenum),function(n) list(n))))
    fraList<-c(fraList,list(sapply((0:prenum),function(n) list(c(1)))))
    {if(k==1)
    {
      return(list(unlist(fraList[[3]])/sum(unlist(fraList[[3]]))))
    }
      else
      {
        for(i in 2:k)
        {
          fraList[[2]]<-lapply(fraList[[1]], function(n) unlist(fraList[[2]][which(fraList[[1]][]==n):(prenum+1)])+n)
          fraList[[3]]<-lapply(fraList[[1]], function(n) unlist(fraList[[3]][which(fraList[[1]][]==n):(prenum+1)]))
          aList<-lapply(fraList[[2]], function(n) unique(n))
          for(j in 1:length(aList))
          {
            fraList[[3]][[j]]<-sapply(aList[[j]], function(n) sum(fraList[[3]][[j]][which(fraList[[2]][[j]][]==n)]))
          }
          fraList[[2]]<-aList
        }
        fraList[[3]]<-sapply(sort(unique(unlist(fraList[[2]]))), function(n) sum(unlist(fraList[[3]])[which(unlist(fraList[[2]])[]==n)]))
        return(list(unlist(fraList[[3]])/sum(unlist(fraList[[3]]))))
      }}
  }
  NumList<-list()
  col<-c()
  conditionClassPre<-conditionList[[1]]
  preNum<-length(conditionList[[1]])
  for(i in 2:length(conditionList))
  {
    col<-c(col,as.vector(sapply(conditionList[[i]],function(n) paste(n,conditionClassPre,sep=","))))
    m=length(conditionList[[i]])
    NumList<-c(NumList,fun_sig(preNum,m))
    preNum<-preNum+m
    conditionClassPre<-c(conditionClassPre,conditionList[[i]])
  }
  ###
  fraV<-NumList[[1]]
  if(length(NumList)>=2)
  {
    for (i in 2:length(NumList)) {
      fraM<-fraV%o%NumList[[i]]
      s<-min(nrow(fraM),ncol(fraM))
      t<-max(nrow(fraM),ncol(fraM))-min(nrow(fraM),ncol(fraM))+1
      fraM<-c(fraM)
      labelV<-c(1:(s-1),rep(s,t),(s-1):1)
      fraV<-c()
      for (j in 1:length(labelV)) {
        fraV<-c(fraV,sum(fraM[1:labelV[j]]))
        fraM<-fraM[-(1:labelV[j])]
      }
    }
  }
  ###
  pr<-c()
  for(i in 1:(floor(length(fraV)/2)))
  {
    pr<-c(pr,fraV[i]+fraV[length(fraV)-i+1])
  }
  if(length(fraV)%%2==1)
  {
    pr<-c(pr,fraV[(length(pr)+1)])
  }
  colAndpr<-list(col,pr)
  return(colAndpr)
}
