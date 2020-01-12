## qdixon.R
# Part of the R package: PMCMRplus
#
# Copyright (C) 2019 Thorsten Pohlert
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#
#' @encoding UTF-8
#' @name Dixon
#' @aliases qdixon
#' @title Dixon distribution
#' @description Density, distribution function, quantile function
#' and random generation for Dixon's ratio statistics \eqn{r_{j,i-1}}{r[j,i-1]}
#' for outlier detection.
#'
#' @details
#' According to McBane (2006) the density of the statistics \eqn{r_{j,i-1}}{r[j,i-1]} of Dixon
#' can be yield if \eqn{x} and \eqn{v} are integrated over the range \eqn{(-\infty < x < \infty, 0 \le v < \infty)}
#'
#' \deqn{
#' \begin{array}[h]{lcl}
#' f(r) & = & \frac{n!}{\left(i-1\right)! \left(n-j-i-1\right)!\left(j-1\right)!} \\
#'       & & \times \int_{-\infty}^{\infty} \int_{0}^{\infty} \left[\int_{-\infty}^{x-v} \phi(t)dt\right]^{i-1}
#'        \left[\int_{x-v}^{x-rv} \phi(t)dt \right]^{n-j-i-1} \\
#'      & &  \times \left[\int_{x-rv}^x \phi(t)dt \right]^{j-1} \phi(x-v)\phi(x-rv)\phi(x)v ~ dv ~ dx \\
#' \end{array}}{%
#'   f(r) = n! / [(i-1)! (n-j-i-1)! (j-1)!]  Int[-\infty,\infty] Int[0,\infty]
#'          [Int[-\infty,x-v] \phi(t)dt ] [Int[x-v,x-rv] \phi(t)dt]^(n-j-i-1)
#'          [Int[x-rv,x] \phi(t)dt]^(j-1) \phi(x-v)\phi(x-rv)\phi(x)v dv dx
#' }
#' where \eqn{v} is the Jacobian and \eqn{\phi(.)} is the density of the standard normal distribution.
#' McBane (2006) has proposed a numerical solution using Gaussian quadratures
#' (Gauss-Hermite quadrature and half-range Hermite quadrature) and coded
#' a library in Fortran. These R functions are wrapper functions to
#' use the respective Fortran code.
#'
#' @note
#' The file \file{slowTest/d-p-q-r-tests.R.out.save} that is included in this package
#' contains some results for the assessment of the numerical accuracy.
#'
#' The slight numerical differences between McBane's original Fortran output
#' (see files \file{slowTests/test[1,2,4].ref.output.txt}) and
#' this implementation are related to different floating point rounding
#' algorithms between R (see \sQuote{round to even} in \code{\link[base]{round}})
#' and Fortran's \code{write(*,'F6.3')} statement.
#'
#'
#' @section Source:
#' The R code is a wrapper to the Fortran code released
#' under GPL >=2 in the electronic supplement of McBane (2006).
#' The original files are \file{rfuncs.f}, \file{utility.f} and \file{dixonr.fi}.
#' They were slightly modified to comply with current CRAN policy and the R manual
#' \sQuote{Writing R Extensions}.
#'
#' @references
#' Dixon, W. J. (1950) Analysis of extreme values.
#' \emph{Ann. Math. Stat.} \bold{21}, 488--506.
#' \url{http://dx.doi.org/10.1214/aoms/1177729747}.
#'
#' Dean, R. B., Dixon, W. J. (1951) Simplified statistics for small
#' numbers of observation. \emph{Anal. Chem.} \bold{23}, 636--638.
#' \url{http://dx.doi.org/10.1021/ac60052a025}.
#'
#' McBane, G. C. (2006) Programs to compute distribution functions
#' and critical values for extreme value ratios for outlier detection.
#' \emph{J. Stat. Soft.} \bold{16}. \url{http://dx.doi.org/10.18637/jss.v016.i03}.
#'
#' @param p vector of probabilities.
# @param n number of observations.
#' @param i number of observations <= x_i
#' @param j number of observations >= x_j
#' @param log.p logical; if \code{TRUE} propabilities p are given as log(p)
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are P[X <= x] otherwise, P[X > x].
#' @return
#' \code{ddixon} gives the density function,
#' \code{pdixon} gives the distribution function,
#' \code{qdixon} gives the quantile function and
#' \code{rdixon} generates random deviates.
#' @keywords distribution
#' @examples
#' set.seed(123)
#' n <- 20
#' Rdixon <- rdixon(n, i = 3, j = 2)
#' Rdixon
#' pdixon(Rdixon, n = n, i = 3, j = 2)
#' ddixon(Rdixon, n = n, i = 3, j = 2)
#'
#' @useDynLib 'dixonTest', .registration = TRUE, .fixes = "F_"
#' @export
qdixon <- function(p, n, i = 1, j = 1, log.p = FALSE, lower.tail = TRUE)
{
    nn <- length(n)
    if (nn > 1) {
        n <- nn
    }

    ok <- checkInput(n, i, j)
    if (!ok) {
        stop("'i','j' must be > 0 and 3 <= n <= 30")
    }

    if(!lower.tail){
        p <- 1 - p
    }
    if (log.p) {
        p <- exp(p)
    }

    m <- as.integer(length(p))
    n <- as.integer(n)
    i <- as.integer(i)
    j <- as.integer(j)
    p <- as.double(p)
    q <- double(m)

    q <-
        .Fortran(
            "forqdixon",
            p = p,
            n = n,
            i = i,
            j = j,
            m = m,
            r = q
        )$r
    return(q)
}

#' @rdname Dixon
#' @aliases pdixon ddixon
#' @param q vector of quantiles
#' @useDynLib 'dixonTest', .registration = TRUE, .fixes = "F_"
#' @keywords distribution
#' @export
pdixon <- function (q, n, i = 1, j = 1, lower.tail = TRUE, log.p = FALSE)
{
    nn <- length(n)
    if (nn > 1) {
        n <- nn
    }

    ok <- checkInput(n, i, j)
    if (!ok) {
        stop("'i','j' must be > 0 and 3 <= n <= 30")
    }

    m <- as.integer(length(q))
    n <- as.integer(n)
    i <- as.integer(i)
    j <- as.integer(j)
    q <- as.double(q)
    p <- double(m)

    p <-
        .Fortran(
            "forpdixon",
            r = q,
            n = n,
            i = i,
            j = j,
            m = m,
            p = p
        )$p

    pval <- sapply(p, function(i)
        min(1, i))
    if (lower.tail) {
        pval <- 1 - pval
    }
    if (log.p) {
        pval <- log(pval)
    }
    return(pval)
}

#' @rdname Dixon
#' @aliases ddixon
#' @param x vector of quantiles.
#' @param log logical; if \code{TRUE} (default),
#' probabilities p are given as log(p).
#' @useDynLib 'dixonTest', .registration = TRUE, .fixes = "F_"
#' @keywords distribution
#' @export
ddixon <- function(x, n, i = 1, j = 1, log = FALSE) {
    nn <- length(n)
    if (nn > 1) {
        n <- nn
    }

    ok <- checkInput(n, i, j)
    if (!ok) {
        stop("'i','j' must be > 0 and 3 <= n <= 30")
    }

    m <- as.integer(length(x))
    n <- as.integer(n)
    i <- as.integer(i)
    j <- as.integer(j)
    x <- as.double(x)
    out <- double(m)

    out <-
        .Fortran(
            "forddixon",
            r = x,
            n = n,
            i = i,
            j = j,
            m = m,
            out = out
        )$out
    if (log) {
        out <- log(out)
    }
    return(out)
}

#' @rdname Dixon
#' @aliases rdixon
#' @param n number of observations. If \code{length(n) > 1},
#' the length is taken to be the number required
#' @importFrom stats runif
#' @export
rdixon <- function(n, i = 1, j = 1){
    nn <- length(n)
    if (nn > 1) {
        n <- nn
    }
    qdixon(runif(n), n, i, j)
}


checkInput <- function(n, i, j) {
    all(n >=3 , n <= 30, i > 0, j > 0, i >= j, i + j < n)
}