#' This script generates draws from a univarite mixture, as described in A. Meucci, "Risk and Asset Allocation",
#' Springer, 2005,  Chapter 4.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 184 - Estimation of a quantile of a mixture I".
#'
#' See Meucci's script for "S_GenerateMixtureSample.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}### 

##################################################################################################################
### Inputs
a   = 0.8;
m_Y = 0.1;
s_Y = 0.2;
m_Z = 0;
s_Z = 0.15;

T   = 52;

##################################################################################################################
### Computations
P = runif(T);
Q = QuantileMixture( P, a, m_Y, s_Y, m_Z, s_Z );

dev.new();
plot( Q );