% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Computing_Row_Probability.R
\name{prow_compute}
\alias{prow_compute}
\title{Perform the step of computing row probability p for all the numbers of 1 under the identified columns for a corresponding submatrix in the ranking matrix}
\usage{
prow_compute(conditionList)
}
\arguments{
\item{conditionList}{the condition partition}
}
\value{
colAndpr: the final identified columns of the corresponding submatrix (its number is denoted by Num(Col)) and all probability p for all the numbers of 1. To meet our destination, we merge them and the result is shown as follows:
the number of 1 (denoted by eNum): 0 and Num(Col), 1 and Num(Col)-1, 2 and Num(Col)-2, ..., Num(Col)/2(or(Num(Col)-1)/2 and (Num(Col)+1)/2)
the probability p:              P(0)+P(Num(Col)),P(1)+P(Num(Col)-1),P(2)+P(Num(Col)-2),...,P(Num(Col)/2)(or P((Num(Col)-1)/2)+P((Num(Col)+1)/2))
only output the probability p,not the number eNum.
}
\description{
Perform the step of computing row probability p for all the numbers of 1 under the identified columns for a corresponding submatrix in the ranking matrix
}
\examples{
data(sim)
aa<-sample(col(sim),3)
bb<-sample(col(sim),2)
conList<-list(aa,bb)
prow_compute(conList)
}
