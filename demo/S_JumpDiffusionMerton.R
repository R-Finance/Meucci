#'This script simulates a jump-diffusion process, as described in A. Meucci, "Risk and Asset Allocation",
#' Springer, 2005,  Chapter 3.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "S_JumoDiffusionMerton.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Parameters
ts = seq( 1/252, 1, 1/252); # grid of time values at which the process is evaluated ("0" will be added, too)
J  = 10; # number of simulations

##################################################################################################################
### Simulate processes

mu  = 0.00; # deterministic drift
sig = 0.20; # Gaussian component

l = 3.45; # Poisson process arrival rate
a = 0;    # drift of log-jump
D = 0.2;  # st.dev of log-jump

X = SimulateJumpDiffusionMerton( mu, sig, l, a, D, ts, J );    
matplot(c( 0, ts), t(X), type="l", xlab = "time", main = "Merton jump-diffusion");
