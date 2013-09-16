\name{MaxRsqTS}
\alias{MaxRsqTS}
\title{Solve for B that maximises sample r-square of F'*B' with X under constraints A*G<=D
and Aeq*G=Deq (A,D, Aeq,Deq conformable matrices),as described in  A. Meucci,
"Risk and Asset Allocation", Springer, 2005.}
\usage{
  MaxRsqTS(X, F, W, A = NULL, D = NULL, Aeq = NULL, Deq,
    lb = NULL, ub = NULL)
}
\arguments{
  \item{X}{: [matrix] (T x N)}

  \item{F}{: [matrix] (T x K)}

  \item{W}{: [matrix] (N x N)}

  \item{A}{: [matrix] linear inequality constraints}

  \item{D}{: [matrix]}

  \item{Aeq}{: [matrix] linear equality constraints}

  \item{Deq}{: [matrix]}

  \item{lb}{: [vector] lower bound}

  \item{ub}{: [vector] upper bound}
}
\value{
  B : [matrix] (N x K)
}
\description{
  Solve for B that maximises sample r-square of F'*B' with
  X under constraints A*G<=D and Aeq*G=Deq (A,D, Aeq,Deq
  conformable matrices),as described in A. Meucci, "Risk
  and Asset Allocation", Springer, 2005.
}
\note{
  Initial code by Tai-Ho Wang
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}. See
  Meucci's script for "MaxRsqTS.m"
}

