#' Power Analysis Tool for Joint Testing Hazards with Competing Risks Data.
#'
#' @param alpha     Type I error.
#' @param beta      Type II error.
#' @param lambda_11 Cause-1 cause-specific hazard in group 1.
#' @param HR_1      Pre-specified cause-1 cause-specific hazard ratio between groups 1 and 2.
#' @param HR_all    Pre-specified any-cause hazard ratio between groups 1 and 2.
#' @param RR        Relative risk of cause-1 failure versus the any-cause failure in group 1.
#' @param r         Length of patient accrual period.
#' @param f         Maximum follow-up period.
#' @param attrition Attrition rate due to lost to follow-up.
#' @param a1        Sample allocation proportion for group 1.
#' @return A dataframe with variables \code{Chi2 Joint}, \code{Maximum Joint}, \code{Bonferroni} methods. The first entry is the required number of cause-1 failures and the second entry is the required total number of patients.
#' @examples
#' library(powerCompRisk)
#' powerCompRisk(alpha = 0.05, beta = 0.2, lambda_11 = 0.3, RR = 0.8,
#' HR_1 = 1.44, HR_all = 1.33, attrition = 0.1, r = 1, f = 8, a1 = 0.5)
#' @references Yang, Q., Fung, W.K., Li, G. (2017) Sample size determination for jointly testing a cause-specific hazard and the any-cause hazard in the presence of competing risks. UCLA Department of Biostatistics Technical Report.
#' @export
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats qnorm pchisq qchisq uniroot
#'
powerCompRisk <- function(alpha, beta, lambda_11, RR, HR_1, HR_all, attrition, r, f, a1) {
  ### Input parameters ###
  alpha      <- alpha
  beta       <- beta
  lambda_11  <- lambda_11
  RR         <- RR
  HR_1       <- HR_1
  HR_all     <- HR_all
  lambda_a1  <- lambda_11 / RR * exp((log(HR_all) - log(HR_1)) / 2)
  lambda_12  <- lambda_11 / HR_1
  lambda_a2  <- lambda_a1 / HR_all
  a2         <- 1 - a1

  # Parameters need to meet a constraint: HR_1 > RR * HR_all
  if (lambda_a2 < lambda_12) {
    print(NA)
    stop("Parameters did not meet the contraint: HR_1 > RR * HR_all")
  }

  # Because attrition=lambda_c/(lambda_c+(lambda_a1+lambda_a2)/2), we have
  lambda_c <- attrition / (1 - attrition) * ((lambda_a1 + lambda_a2) / 2)

  ##########################################################################
  ## Number of cause-1 failures required by the maximum joint test
  ##########################################################################
  ## Find c_m by solving the alpha equation
  corr <- matrix(c(1, sqrt(RR), sqrt(RR), 1), nrow = 2, ncol = 2)

  cv <- function(x) {
    (1 - pmvnorm(c(-x, -x), c(x, x), c(0, 0), corr)) - alpha
  }

  C_m <- uniroot(cv, c(0, 3), tol = 0.0001)$root

  ## Find D by solving the beta/power equaiton
  cv2 <-function(x) {
    beta - pmvnorm(c(-C_m, -C_m), c(C_m, C_m), c(log(HR_1) * sqrt(a1 * a2 * x),
                                               log(HR_all)  *sqrt(a1 * a2 * x / RR)), corr)
  }
  d_max <- ceiling(uniroot(cv2, c(0, 5000), tol = 0.0001)$root)

  ####################################################
  ## Number of cause-1 failures required by the chi-square joint test
  ####################################################
  ## Find D, number of event of interest, for the chisquare joint test by iterations

  C <- qchisq(1 - alpha, 2, 0)

  cv_chisq <- function(x) {
    beta - pchisq(C, 2, ncp = x * a1 * a2 * (log(HR_1)^2 - 2 * log(HR_1) * log(HR_all) + log(HR_all)^2 / RR) / (1 - RR))
  }

  d_chi <- ceiling(uniroot(cv_chisq, c(0, 5000), tol = 0.0001)$root)

  ##################################################
  ## Number of cause-1 failures required by the Bonferroni method
  ##################################################
  # CSH_1
  d1 <- ceiling((qnorm(1 - alpha / 4) + qnorm(1 - beta))^2 / (a1 * a2 * log(HR_1)^2) )
  # ACH
  d2 <- ceiling((qnorm(1 - alpha / 4) + qnorm(1 - beta))^2 / (a1 * a2 * log(HR_all)^2))

  lambda_all1 <- lambda_a1 + lambda_c
  lambda_all2 <- lambda_a2 + lambda_c

  # Calculate the probability of observing a cause-1 failure by the end of study
  if (r == 0) {
    P_11 <- lambda_11 / lambda_all1 * (1 - (exp(-lambda_all1 * f)))
    P_12 <- lambda_12 / lambda_all2 * (1 - (exp(-lambda_all2 * f)))
    P_a1 <- lambda_a1 / lambda_all1 * (1 - (exp(-lambda_all1 * f)))
    P_a2 <- lambda_a2 / lambda_all2 * (1 - (exp(-lambda_all2 * f)))
  } else {
    P_11 <- lambda_11 / lambda_all1 * (1 - (exp(-lambda_all1 *(f -r)) - exp(-lambda_all1 * f)) / lambda_all1 * r)
    P_12 <- lambda_12 / lambda_all2 * (1 - (exp(-lambda_all2 *(f -r)) - exp(-lambda_all2 * f)) / lambda_all2 * r)
    P_a1 <- lambda_a1 / lambda_all1 * (1 - (exp(-lambda_all1 *(f -r)) - exp(-lambda_all1 * f)) / lambda_all1 * r)
    P_a2 <- lambda_a2 / lambda_all2 * (1 - (exp(-lambda_all2 *(f -r)) - exp(-lambda_all2 * f)) / lambda_all2 * r)
  }

  P_1 <- a1 * P_11 + a2 * P_12
  P_a <- a1 * P_a1 + a2 * P_a2
  n1  <- ceiling(d1 / P_1)
  n2  <- ceiling(d2 / P_a)
  ##############################################################################
  # Calculate the required number of patietns for each of the three joint tests
  ##############################################################################
  # Required sample size for the chi-square joint test
  samplesize_chi <- ceiling(d_chi / P_1)
  # Required sample size for the maximum joint test
  samplesize_bon <- ceiling(min(n1, n2))
  # Required sample size for the Bonferroni joint test
  samplesize_max <- ceiling(d_max / P_1)
  d_bon <- ceiling(min(n1, n2) * P_1)

  #################################
  # Output results
  #################################
  result <- matrix(nrow = 2, ncol = 3)
  result[1, 1] <- ceiling(d_chi / 2) * 2
  result[2, 1] <- ceiling(samplesize_chi / 2) * 2
  result[1, 2] <- ceiling(d_max / 2) * 2
  result[2, 2] <- ceiling(samplesize_max / 2) * 2
  result[1, 3] <- ceiling(d_bon / 2) * 2
  result[2, 3] <- ceiling(samplesize_bon / 2) * 2

  writeLines(paste0("Number of death for Chi-square Joint Test: ", d_chi))
  writeLines(paste0("Total number of patients for Chi-square Joint Test: ", samplesize_chi, "\n"))
  writeLines(paste0("Number of death for Maximum Joint Test: ", d_max))
  writeLines(paste0("Total number of patients for Maximum Joint Test: ", samplesize_max, "\n"))
  writeLines(paste0("Number of death for Bonferroni Test: ", d_bon))
  writeLines(paste0("Total number of patients for Bonferroni Test: ", samplesize_bon, "\n"))

  result <- as.data.frame(result)
  colnames(result) <- c('Chi2 Joint', 'Maximum Joint', 'Bonferroni')
  rownames(result) <- c('D1', 'N')
  return(result)
}
