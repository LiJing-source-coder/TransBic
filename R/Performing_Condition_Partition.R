#'
#'Perform the step of performing condition partition for a set of filtered ordered pairs of condition
#'
#'@param conditionAll all conditions of the input expression matrix
#'@param Rfilt the filtered R, which R is the primarily identified columns of a corresponding submatrix in the ranking matrix
#'@param mfra the predefined parameter, and refer to "TransBic.R"
#'@param cfra the predefined parameter, and refer to "TransBic.R"
#'@param minC the predefined parameter, and refer to "TransBic.R"
#'
#'@return conditionList: the resulting condition partition according to the input Rfilt.
#'
#'@export
#'
#'@examples
#'data(sim)
#'aa<-sample(col(sim),3)
#'bb<-sample(col(sim),2)
#'Rex<-as.vector(sapply(aa,function(x) paste(x,bb,sep =",")))
#'condition_partition(col(sim),Rex,mfra=0.8,cfra=0.75,minC=2)
#'
"condition_partition"<-function(conditionAll,Rfilt,mfra,cfra,minC)
{
  comSet_compute<-function(x,comSet)
  {
    coml<-length(intersect(names(which(x==1)),comSet))
    return(coml)
  }
  ###
  if(length(Rfilt)>=1)
  {adj_matrix<-matrix(0,nrow = length(conditionAll),ncol = length(conditionAll),dimnames = list(conditionAll,conditionAll))
  for (i in 1:length(Rfilt)) {
    conPair<-unlist(strsplit(Rfilt[i],split = ","))
    adj_matrix[conPair[1],conPair[2]]<-1
  }
  ###
  canNod<-colnames(adj_matrix)
  conditionList<-list()
  conditionRe<-colnames(adj_matrix)
  while(length(canNod)>0)
  {
    keyNod<-canNod[which.max(apply(as.matrix(adj_matrix[conditionRe,canNod]),2,sum))]
    botcomSet<-names(which(adj_matrix[,keyNod]==1))
    topcomSet<-names(which(adj_matrix[keyNod,]==1))
    ###
    botInSet<-apply(as.matrix(adj_matrix[,conditionRe]),2,comSet_compute,comSet=botcomSet)
    topInSet<-apply(as.matrix(adj_matrix[conditionRe,]),1,comSet_compute,comSet=topcomSet)
    ###
    concluster<-c()
    for(i in 1:length(conditionRe))
    {
      afra<-(as.double(topInSet[i])+as.double(botInSet[i]))/(length(topcomSet)+length(botcomSet))
      if(afra>mfra)
      {
        concluster<-c(concluster,conditionRe[i])
      }
    }
    conditionList<-c(conditionList,list(concluster))
 ###
    topCluster<-unlist(conditionList[1:length(conditionList)])
    conditionRe<-setdiff(colnames(adj_matrix),topCluster)
    if(length(conditionRe)<=1)
    {break}
    canNod<-conditionRe[which(apply(as.matrix(adj_matrix[conditionRe,topCluster]),1,sum)>cfra*length(topCluster))]
  }
 ###
  clusterNum<-sapply(conditionList,function(n) length(n))
  while(length(conditionList)>1 && (minC<length(conditionAll)) && (length(which(clusterNum<minC))>0))
  {
    merC<-which(clusterNum<minC)
    merRat<-c()
    parC<-c()
    for(i in merC)
    {
      headSet<-setdiff(unlist(conditionList[1:i]),unlist(conditionList[[i]]))
      tailSet<-setdiff(unlist(conditionList[i:length(conditionList)]),unlist(conditionList[[i]]))
      headInSet<-sum(adj_matrix[unlist(conditionList[[i]]),headSet])
      tailInSet<-sum(adj_matrix[tailSet,unlist(conditionList[[i]])])
      {if(length(headSet)==0)
      {headRat<-0
      parC<-c(parC,"t")
      merRat<-c(merRat,tailInSet/(length(tailSet)*length(unlist(conditionList[[i]]))))}
      else if(length(tailSet)==0)
      {tailRat<-0
      parC<-c(parC,"h")
      merRat<-c(merRat,headInSet/(length(headSet)*length(unlist(conditionList[[i]]))))}
        else
        {headRat<-headInSet/(length(headSet)*length(unlist(conditionList[[i]])))
        tailRat<-tailInSet/(length(tailSet)*length(unlist(conditionList[[i]])))
        htRat<-c(headRat,tailRat)
        names(htRat)<-c("h","t")
        merRat<-c(merRat,min(htRat))
        parC<-c(parC,names(which.min(htRat)))}}}
    if(parC[which.min(merRat)]=="h")
    {conditionList[merC[which.min(merRat)]]<-list(c(unlist(conditionList[[merC[which.min(merRat)]]]),unlist(conditionList[[merC[which.min(merRat)]-1]])))
    conditionList<-conditionList[-(merC[which.min(merRat)]-1)]}
    if(parC[which.min(merRat)]=="t")
    {conditionList[merC[which.min(merRat)]]<-list(c(unlist(conditionList[[merC[which.min(merRat)]]]),unlist(conditionList[[merC[which.min(merRat)]+1]])))
    conditionList<-conditionList[-(merC[which.min(merRat)]+1)]}
    clusterNum<-sapply(conditionList,function(n) length(n))
    }}
  return(conditionList)
  }



