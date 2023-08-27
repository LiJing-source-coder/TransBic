#'
#'Perform the step of constructing ranking matrix for the input gene expression matrix
#'
#'@param data the input gene expression matrix
#'       (each row represents a gene and each column represents a condition)
#'@return M_matrix: the corresponding ranking matrix
#'
#'@export
#'
#'@examples
#'data(sim)
#'rankingM_constructing(sim)
#'

"rankingM_constructing"<-function(data)
{
  matrix_new<-function(x,c_matrix)
  {
    a<-x[1]
    b<-x[2]
    l<-nrow(c_matrix)
    newVec1<-t(t(rep(0,2*nrow(c_matrix))))
    rowNum<-which(as.numeric(c_matrix[,a])<as.numeric(c_matrix[,b]))
    rowsetNum<-which(as.numeric(c_matrix[,a])>as.numeric(c_matrix[,b]))
    if(length(rowNum)!=0)
    {newVec1[rowNum,]<-1}
    if(length(rowsetNum)!=0)
    {newVec1[l+rowsetNum,]<-1}
    return(newVec1)
  }
  labelCol<-t(utils::combn(colnames(data),2))
  M_matrix<-matrix(0,nrow = 2*nrow(data),ncol = 2*nrow(labelCol))
  M_matrix[,1:nrow(labelCol)]<-as.matrix(apply(labelCol,1,matrix_new,c_matrix=data))
  M_matrix[,(nrow(labelCol)+1):ncol(M_matrix)]<-as.matrix(rbind(M_matrix[(nrow(data)+1):nrow(M_matrix),1:nrow(labelCol)],M_matrix[1:nrow(data),1:nrow(labelCol)]))
  labelCol<-rbind(labelCol,t(apply(labelCol,1,rev)))
  labelColNam<-apply(labelCol,1,function(x) paste(x,collapse =","))
  colnames(M_matrix)<-labelColNam
  labelRowNam<-paste0(row.names(data),"_1")
  labelRowNam<-c(labelRowNam,paste0(row.names(data),"_2"))
  rownames(M_matrix)<-labelRowNam
  M_matrix<-as.matrix(M_matrix)
  return(M_matrix)
}
