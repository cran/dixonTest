# dixonTest.R
# Part of the R package: dixonTest
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

#' @name dixonTest
#' @title Dixons Outlier Test (Q-Test)
#' @description
#' Performs Dixons single outlier test.
#'
#' @details
#' Let \eqn{X} denote an identically and independently distributed
#' normal variate. Further, let the increasingly ordered realizations
#' denote \eqn{x_1 \le x_2 \le \ldots \le x_n}.
#' Dixon (1950) proposed the following ratio statistic to detect
#' an outlier (two sided):
#'
#'  \deqn{
#'   r_{j,i-1} = \max\left\{\frac{x_n - x_{n-j}}{x_n - x_i},
#'   \frac{x_{1+j} - x_1}{x_{n-i} - x_1}\right\}}{%
#'   r[j,i-1] = max\{(x[n] - x[n-j]) / (x[n] - x[i]),
#'              (x[1+j] - x[1]) / (x[n-i] - x[1])\}
#'   }
#'
#' The null hypothesis, no outlier, is tested against the alternative,
#' at least one observation is an outlier (two sided). The subscript \eqn{j}
#' on the \eqn{r} symbol indicates the number of
#' outliers that are suspected at the upper end of the data set,
#' and the subscript \eqn{i} indicates the number of outliers suspected
#' at the lower end. For \eqn{r_{10}} it is also common to use the
#' statistic \eqn{Q}.
#'
#' The statistic for a single maximum outlier is:
#' \deqn{
#'   r_{j,i-1} = \left(x_n - x_{n-j} \right) / \left(x_n - x_i\right)}{%
#'   r[j,i-1] = (x[n] - x[n-j]) / (x[n] - x[i])
#'   }
#' The null hypothesis is tested against the alternative,
#' the maximum observation is an outlier.
#'
#'
#' For testing a single minimum outlier, the test statistic is:
#' \deqn{
#'    r_{j,i-1} = \left(x_{1+j} - x_1 \right) / \left(x_{n-i} - x_1 \right)}{%
#'    r[j,i-1] = (x[1+j] - x[1]) / (x[n] - x[i])
#'    }
#'
#' The null hypothesis is tested against the alternative,
#' the minimum observation is an outlier.
#'
#' Apart from the earlier Dixons Q-test (i.e. \eqn{r_{10}}),
#' a refined version that was later proposed by Dixon can be performed
#' with this function, where the statistic \eqn{r_{j,i-1}} depends on
#' the sample size as follows:
#'
#' \tabular{rl}{
#'  \eqn{r_{10}}: \tab \eqn{3 \le n \le 7} \cr
#'  \eqn{r_{11}}: \tab \eqn{8 \le n \le 10} \cr
#'  \eqn{r_{21}}; \tab \eqn{11 \le n \le 13} \cr
#'  \eqn{r_{22}}: \tab \eqn{14 \le n \le 30} \cr
#' }
#'
#' The p-value is computed with the function \code{\link{pdixon}}.
#' @inherit Dixon references
#' @param x a numeric vector of data
#' @param alternative the alternative hypothesis.
#' Defaults to \code{"two.sided"}
#' @param refined logical indicator, whether the refined version
#' or the Q-test shall be performed. Defaults to \code{FALSE}
#' @keywords htest
#' @concept outliers
#' @examples
#' ## example from Dean and Dixon 1951, Anal. Chem., 23, 636-639.
#' x <- c(40.02, 40.12, 40.16, 40.18, 40.18, 40.20)
#' dixonTest(x, alternative = "two.sided")
#'
#' ## example from the dataplot manual of NIST
#' x <- c(568, 570, 570, 570, 572, 578, 584, 596)
#' dixonTest(x, alternative = "greater", refined = TRUE)
#'
#' @importFrom stats na.omit
#' @export
dixonTest <- function(x, alternative = c("two.sided", "greater", "less"), refined = FALSE){

    if (!is.numeric(x)){
        stop("'x' must be a numeric vector.")
    }
    alternative <- match.arg(alternative)
    data.name <- deparse(substitute(x))
    x <- na.omit(x)

    n <- length(x)
    if (n < 3 | n > 30) {
        stop("'x' must be a vector of length 2 < n 30")
    }

    o <- order(x)
    xx <- x[o]

    ## generic function
    if (!refined | 3 <=n & n <= 7) {
        i <- 1
        j <- 1
    } else if (8 <= n & n <= 10) {
        i <- 2
        j <- 1
    } else if (11 <= n & n <= 13) {
        i <- 2
        j <- 2
    } else if (14 <= n & n <= 30) {
        i <- 3
        j <- 2
    }

    up <- function(n, i, j) {
        (xx[n] - xx[n-j]) /
            (xx[n] - xx[i])
    }

    lo <- function(n, i, j) {
        (xx[1+j] - xx[1]) /
            (xx[n] - xx[i])
    }

    if (alternative == "two.sided") {
    ## maximumor minimum is outlier?
        u <- up(n,i,j)
        l <- lo(n,i,j)
        tmp <- c(u,l)
        w <- which.max(tmp)
        ii <- ifelse(w == 1, which.max(x), which.min(x))
        rji <- tmp[w]

    } else if (alternative == "less") {
        ## testing for outlying minimum
        ii <- which.min(x)
        rji <- lo(n, i, j)

    } else {
        ## testing for outlying maximum
        ii <- which.max(x)
        rji <- up(n, i, j)

    }
    i <- as.integer(i)
    j <- as.integer(j)
    n <- as.integer(n)
    rji <- as.double(rji)


    ## calculate p-value
    pval <- pdixon(rji, n, i, j, lower.tail = TRUE)
    if (alternative == "two.sided") {
        pval <- min(1, 2 * pval)
    }

    val <- x[ii]
    Stat <- rji
    names(Stat) <- ifelse(!refined, "Q", paste0("r",j, i-1))

    ans <- list(
        method = "Dixons outlier test",
        alternative = alternative,
        statistic = Stat,
        parameter = c("n" = n),
        p.value = pval,
        estimate = c(c("pos" = ii),
                     c("value" = val)),
        data.name = data.name)
    class(ans) <- "htest"
    return(ans)
}
