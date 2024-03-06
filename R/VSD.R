#' Hi! Welcome to VSD, a tool for calculating viral replicons and viral particles abundance
#' based on viral RNA expression profiles.
#' @title  Viral Shedding Detector.
#'
#' @description  VSD is an optimization model developed using a linear programming algorithm which includes an objective function, constraints, and decision variables. VSD, considers the abundance of viral replicons and viral particles as variables in its optimization process.
#'
#' @details  Input data:  Normalized viral RNA sequencing data, and the normalization should be corrected for RNA length;
#' Output data: Abundance of viral replicons, viral particles and viral fragments (fitted residual representation)
#'
#'@param expr A matrix or data frame consisting of viral RNA expression values.
#'@param reference Simulated viral RNA detection patterns of viruses inside or outside the cell.
#'
#'@usage VSD(expr, reference)
#'
#'@aliases VSD
#'
#' @return A dataframe which contains three parts: viral replicons, viral particles and viral fragments abundance.
#'
#'@export VSD
#'@examples VSD(expr=example,reference=ref)
#
#'@importFrom Rglpk Rglpk_solve_LP
#
#
#
#
#
VSD<-function(expr,reference){
  jieguo_lpm<-c()
  for (i in 1:dim(expr)[2]) {

    obj<-c(rep(1,dim(reference)[2]))
    mat<-reference
    rhs<-expr[,i]
    dir<-c(rep("<=",dim(reference)[1]))
    max<-T
    if(requireNamespace("Rglpk",quietly = T)){
    a1<-Rglpk::Rglpk_solve_LP(obj = obj,mat = mat,dir = dir,rhs = rhs,max = max)$solution
    }
    jieguo_lpm<-cbind(jieguo_lpm,a1)
  }
  rownames(jieguo_lpm)<-colnames(reference)
  colnames(jieguo_lpm)<-colnames(expr)
  nihejieguo<-as.matrix(reference)%*%jieguo_lpm
  cancha<-expr-nihejieguo
  residuals<-apply(abs(cancha),2,sum)
  jieguo_lpm<-as.data.frame(rbind(jieguo_lpm,residuals))
  return(jieguo_lpm)
}
