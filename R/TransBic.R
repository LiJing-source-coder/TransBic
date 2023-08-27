#'
#'Perform the TransBic program
#'
#'@param data normalized gene expression matrix(each column represents a condition
#'and each row represents a gene)
#'@param minC minimum number of conditions in each condition cluster from resulting condition partitions(int)
#'@param minR minimum number of rows of resulting biclusters(int)
#'@param e0 error rate of the identified rows in the corresponding submatrix in the ranking matrix([0,1])
#'@param tfra consistency level of the identified columns of the corresponding submatrix in this iteration([0,1])
#'compared to the last iteration([0,1])
#'@param tNum number of consecutive iterations for keeping the above consistence(int)
#'@param cfra level to evaluate whether a selected condition is from a new condition cluster from the same partition with the
#'identified condition clusters ([0,1])
#'@param mfra level to evaluate whether two conditions are from the same condition cluster([0,1])
#'@param sfra fraction of the identified columns to preserve as the columns of the corresponding submatrix([0,1])
#'
#'@return Biclusters
#'Taken a bicluster as an example, its output format is shown as follows:
#'line 1_rank order of biclusters:        Bic:number
#'line 2_rows:                            genes_1/2
#'line 3_columns:                         condition partition for a identified condition subset
#'        (Suppose the condition partition is C1,C2,C3,...,CT from the top down,and then
#'        the row g_1 has [g,C1]> [g,C2]>[g,C3]>⋯>[g,CT], where [g,C] denote all expression values under condition cluster C for gene g.
#'        At the same time,  the row g_2 has [g,C1]< [g,C2]<[g,C3]<⋯<[g,CT])
#'
#'@export
#'
#'@examples
#'data(sim)
#'TransBic(sim,tNum=3,sfra = 0.9)
TransBic<-function(data,minC=2,minR=2,e0=0.1,tfra=0.75,tNum=5,cfra=0.75,mfra=0.8,sfra=1)
{
  fraV_compute=function(c)
{
  fraV<-sum(c)/length(c)
  return(fraV)
  }
M_matrix<-rankingM_constructing(data)
K_matrix<-M_matrix%*%t(M_matrix)
bicList<-list()
bicNum=1
sNum<-0
conditionAll<-colnames(data)
while ((nrow(K_matrix)>1) && (max(K_matrix[lower.tri(K_matrix)])>0)){
  maxFre<-max(K_matrix[lower.tri(K_matrix)])
  seed<-which(K_matrix==maxFre,arr.ind = T)[1,]
  bicDf<-data.frame(bicRow=row.names(K_matrix)[seed[1]],canRow=colnames(K_matrix)[seed[2]],canNum=maxFre,stringsAsFactors=FALSE)
  comRow<-bicDf$canRow
  rowdex<-bicDf$bicRow
  Num<-0
  Del<-0
  bicExi<-0
  pNum<-0
  R_last<-c()
  repeat{
    bicDf<-rbind(bicDf,data.frame(bicRow=comRow,canRow="0",canNum=0))
    rowdex<-union(rowdex,comRow)
    if(length(rowdex)==nrow(K_matrix))
    {break}
    canS<-c(bicDf$bicRow[which(bicDf$canRow==comRow)],comRow)
    cc_matrix<-K_matrix[canS,setdiff(colnames(K_matrix),rowdex), drop=FALSE]
    cNum<-apply(as.matrix(cc_matrix),1,max)
    bicDf[which(bicDf$bicRow%in%canS),"canNum"]<-as.integer(cNum)
    cRow<-apply(cc_matrix,1,which.max)
    bicDf[which(bicDf$bicRow%in%canS),"canRow"]<-as.character(colnames(cc_matrix)[cRow])
    {if(Del==1)
    {
      bicDf<-bicDf[-nrow(bicDf),]
      Num<-Num+1
    }
      else
      {
        Num<-0
      }}
    if(Num>20)
    {break}
    if(Del==0)
    {
      fraV<-apply(as.matrix(M_matrix[bicDf$bicRow,]),2,fraV_compute)
      flag<-try(stats::kmeans(fraV,3),silent=TRUE)
      {if('try-error' %in% class(flag))
      {km<-stats::kmeans(fraV,2)
      colClass<-list()
      for (m in 1:2) {
        colClass<-c(colClass,list(names(which(km[["cluster"]]==m))))
      }
      colClassNum<-order(c(mean(fraV[colClass[[1]]]),mean(fraV[colClass[[2]]]),decreasing = T))
      R<-colClass[[colClassNum[1]]]
      Q<-colClass[[colClassNum[2]]]}
      else
      {km<-stats::kmeans(fraV,3)
      colClass<-list()
      for (m in 1:3) {
        colClass<-c(colClass,list(names(which(km[["cluster"]]==m))))
      }
      colClassNum<-order(c(mean(fraV[colClass[[1]]]),mean(fraV[colClass[[2]]]),mean(fraV[colClass[[3]]])),decreasing = T)
      R<-colClass[[colClassNum[1]]]
      Q<-colClass[[colClassNum[3]]]}}
      if(length(R_last)==0)
      {R_last<-R}
      {if((length(intersect(R,R_last))/max(length(R),length(R_last)))>tfra)
      {
        pNum<-pNum+1
      }
        else
        {
          pNum<-0
        }}
      R_last<-R
    }
    if(pNum>tNum)
    {break}
    addRow<-bicDf[which.max(bicDf$canNum),"canRow"]
    p1<-length(which(M_matrix[addRow,R]==1))
    q1<-length(which(M_matrix[addRow,R]==0))
    p2<-length(which(M_matrix[addRow,Q]==1))
    q2<-length(which(M_matrix[addRow,Q]==0))
    {if((p1==0) || (q2==0))
    {Del<-1}
      else if((q1==0) || (p2==0))
      {Del<-0}
      else
      {kl<-p1/(p1+q1)*log((p1/(p1+q1))/(p2/(p2+q2)))+q1/(p1+q1)*log((q1/(p1+q1))/(q2/(p2+q2)))
      {if(kl<0.3)
      {Del<-1}
        else
        {Del<-0}}}}
    comRow<-addRow
  }
  ####
  if(length(rowdex)==nrow(K_matrix))
  {
    if(Del==0)
    {fraV<-apply(as.matrix(M_matrix[bicDf$bicRow,]),2,fraV_compute)
    flag<-try(stats::kmeans(fraV,3))
    {if('try-error' %in% class(flag))
    {km<-stats::kmeans(fraV,2)
    colClass<-list()
    for (m in 1:2) {
      colClass<-c(colClass,list(names(which(km[["cluster"]]==m))))
    }
    colClassNum<-order(c(mean(fraV[colClass[[1]]]),mean(fraV[colClass[[2]]]),decreasing = T))
    R<-colClass[[colClassNum[1]]]
    Q<-colClass[[colClassNum[2]]]}
      else
      {km<-stats::kmeans(fraV,3)
      colClass<-list()
      for (m in 1:3) {
        colClass<-c(colClass,list(names(which(km[["cluster"]]==m))))
      }
      colClassNum<-order(c(mean(fraV[colClass[[1]]]),mean(fraV[colClass[[2]]]),mean(fraV[colClass[[3]]])),decreasing = T)
      R<-colClass[[colClassNum[1]]]
      Q<-colClass[[colClassNum[3]]]}}
    {if((length(intersect(R,R_last))/max(length(R),length(R_last)))>tfra)
    {
      pNum<-pNum+1
    }
      else
      {
        pNum<-0
      }}
    }
    else
    {
      bicDf<-bicDf[-nrow(bicDf),]
    }
  }
#######
 comRowdel<-bicDf$bicRow
 if(pNum>tNum)
  {
    scol2<-M_matrix[which(row.names(M_matrix)%in%bicDf[,1]),which(colnames(M_matrix)%in%R)]
    dfra<-apply(scol2, 2, sum)
    dfra<-dfra/nrow(scol2)
    Rorder<-order(dfra,decreasing=TRUE)[1:floor(length(dfra)*sfra)]
    Rfilt<-names(dfra)[Rorder]
    conditionList<-condition_partition(conditionAll,Rfilt,mfra,cfra,minC)
    if(length(conditionList)>=2)
    {
      res<-prow_compute(conditionList)
      resCol<-res[[1]]
      resPr<-res[[2]]
      rowfra<-ceiling(length(resCol)*(1-e0))
      rowNum<-apply(M_matrix[,resCol],1,sum)
      p<-sum(resPr[1:(length(resCol)-rowfra+1)])
      rowNa<-names(which(rowNum>=rowfra))
            if(length(rowNa)>=minR)
      {s=1-stats::pbinom(length(rowNa)-1,nrow(M_matrix)/2,p)
      if(s<0.001)
      {bicList<-c(bicList,list(paste("Bic:",bicNum,sep = "")))
      bicList<-c(bicList,list(rowNa))
      bicList<-c(bicList,list(conditionList))
      bicNum=bicNum+1
      bicExi<-1
      }}
      comRowdel<-rowNa
    }
 }
  {if(bicExi==0)
  {
    sNum<-sNum+1
  }
    else
    {
      sNum<-0
    }}
  if(sNum>15)
  {break}
    if(bicExi==1)
  {comRowdel<-substr(comRowdel,1,nchar(comRowdel)-2)
  rowDel<-which(substr(row.names(K_matrix),1,nchar(row.names(K_matrix))-2)%in%comRowdel)
  K_matrix<-K_matrix[-rowDel,-rowDel]}
  if(bicExi==0)
  {K_matrix[row.names(K_matrix)%in%comRowdel,colnames(K_matrix)%in%comRowdel]<-0
  rowEndf<-substr(comRowdel,nchar(comRowdel)-2,nchar(comRowdel))
  rowEnd<-sapply(rowEndf, function(n) setdiff(c("_1","_2"),n))
  rowCle<-paste0(substr(comRowdel,1,nchar(comRowdel)-2),rowEnd)
  K_matrix[row.names(K_matrix)%in%rowCle,colnames(K_matrix)%in%rowCle]<-0}
}

{if(length(bicList)==0)
  print("No Bicluster!")
  else
  {if(length(bicList)>3)
  {
    i=1
    while(i<(length(bicList)/3))
    {
      bicdel<-c()
      for(j in (i+1):(length(bicList)/3))
      {
        if(length(setdiff(bicList[[(3*i-1)]],bicList[[(3*j-1)]]))==0)
        {bicdel<-c(bicdel,j)}
      }
      bicDel<-c(3*bicdel-2,3*bicdel-1,3*bicdel)
      if(length(bicDel)>0)
      {bicList<-bicList[-bicDel]}
      i=i+1
    }
  }
    for(n in 1:(length(bicList)/3))
    {
      bicList[[3*n-2]]<-paste("Bic:",n,sep = "")
    }
  }}
return(bicList)}

#'a simulated dataset
#'
#'In the dataset, a bicluster with rows {1,2,3,4,5,6},columns {1,2,3,4,5,6} and
#'its pattern: for g={1,3,5,6}, [g,{1,2}]>[g,{3,4}]>[g,{5,6}] and for g={2,4}, [g,{1,2}]<[g,{3,4}]<[g,{5,6}]is implanted.
#'
#'@format A matrix with 8 rows and 8 columns:
#'
#'@usage data(sim)
#'
#'@examples
#'#data(sim)
#'#View(sim)
"sim"
