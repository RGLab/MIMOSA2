#' Test for a difference of odds ratios
#' @details
#' SE of log odds = sqrt(1/n1+1/n2+1/n3+1/n4)
#' SE of difference in log odds  = sqrt(SE1^2 +SE2^2)
#' delta  = difference in log odds
#' Test delta / SE from standard normal.
#'
#' @param Ntot \code{matrix} integer vector of total trials. One row per subject. Should have four columns named "ns1" "ns0" "nu1" and "nu0"
#' @param ns1 \code{numeric} integer vector of successes in condition 1 treatment s.
#' @param nu1 \code{numeric} integer vector of successes in condition 1 treatment u.
#' @param ns0 \code{numeric} integer vector of successes in condition 0 treatment s.
#' @param nu0 \code{numeric} integer vector of successes in condition 0 treatment u.
#' @return \code{numeric} vector of one-sided p-values from a standard normal.
#' @export
#' @seealso \link{MIMOSA2}
#' @examples
#' s = simulate_MIMOSA2()
#' R = MIMOSA2(Ntot=s$Ntot, ns1 = s$ns1, nu1 = s$nu1, nu0 = s$nu0, ns0 = s$ns0)
#' or_result = ORTest(Ntot = s$Ntot, ns1 = s$ns1, nu1 = s$nu1, nu0 = s$nu0, ns0 = s$ns0)
ORTest = function(Ntot, ns1, nu1, ns0, nu0) {
  or_test = sapply(1:nrow(Ntot), function(i) {
    or_delta = log(fisher.test(rbind(
      cbind(Ntot[i, "ns1"] + 1, ns1[i] + 1),
      cbind(Ntot[i, "nu1"] + 1, nu1[i] + 1)
    ))$estimate) -   log(fisher.test(rbind(
      cbind(Ntot[i, "ns0"] + 1, ns0[i] + 1),
      cbind(Ntot[i, "nu0"] +
              1, nu0[i] + 1)
    ))$estimate)
    se_or_delta =  sqrt(sum(sqrt(sum(
      1 / c(cbind(Ntot[i, "ns1"] + 1, ns1[i] + 1),
            cbind(Ntot[i, "nu1"] + 1, nu1[i] + 1))
    )) ^ 2,
    sqrt(sum(
      1 / c(cbind(Ntot[i, "ns0"] + 1, ns0[i] + 1),
            cbind(Ntot[i, "nu0"] + 1, nu0[i] + 1))
    )) ^ 2))
    pnorm(or_delta / se_or_delta) #one sided test
  })
  return(or_test)
}

#'Generate matrix for ROC plot.
#'
#' @param or_test output of \code{ORTest}
#' @param z output of MIMOSA2 \code{z} component
#' @param truth \code{vector} of ground truth of length number of observations where response = "R".
#' @return \code{matrix} with columns "FPR" "TPR" and "Method"
#' @export
#' @seealso \link{MIMOSA2}
#' @examples 
#' s = simulate_MIMOSA2()
#' R = MIMOSA2(Ntot=s$Ntot, ns1 = s$ns1, nu1 = s$nu1, nu0 = s$nu0, ns0 = s$ns0)
#' or_result = ORTest(Ntot = s$Ntot, ns1 = s$ns1, nu1 = s$nu1, nu0 = s$nu0, ns0 = s$ns0)
#' ROC(or_test = or_result, z = R$z, truth = s$truth)
ROC = function(or_test, z, truth) {
  toplot = data.table(rbind(
    roc(or_test, truth == "R"),
    roc(1 - rowSums(z[, c(1, 5)]), truth == "R")
  ),
  Method = gl(
    n = 2,
    k = 1000,
    labels = c("Test of Difference in Log-odds Ratio",
               "MIMOSA with baseline")
  ))
  return(toplot)
}


#' Plot ROCs using ggplot
#' @name ROCPlot
#' @details uses the output of ROC
#' @param R \code{matrix} output of \code{ROC()}
#' @export
#' @seealso \link{MIMOSA2}
#' @examples 
#' s = simulate_MIMOSA2()
#' R = MIMOSA2(Ntot=s$Ntot, ns1 = s$ns1, nu1 = s$nu1, nu0 = s$nu0, ns0 = s$ns0)
#' or_result = ORTest(Ntot = s$Ntot, ns1 = s$ns1, nu1 = s$nu1, nu0 = s$nu0, ns0 = s$ns0)
#' roc = ROC(or_test = or_result, z = R$z, truth = s$truth)
#' ROCPlot(roc)
ROCPlot = function(R) {
  ggplot(R) +
    aes_string(x = "FPR", y = "TPR", color = "Method") +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_color_brewer(type = "qual", palette = 6) +
    scale_x_continuous("False Positive Rate") +
    scale_y_continuous("True Positive Rate") +
    stat_quantile(
      method = "rqss",
      quantiles = c(0.05, 0.95),
      linetype = "dashed",
      lambda = 1
    ) + stat_quantile(method = "rqss",
                      quantiles = c(0.5),
                      lambda = 1) + theme(axis.text.x = element_text(angle = 45))
}

#' Boxplots from MIMOSA fit
#'
#' @param obj \code{MIMOSA2()} output
#' @param truth \code{character vector} of true responder / non-responder calls. Responder should be "R"
#' @return ggplot grob
#' @export
#' @seealso \link{MIMOSA2}
#' @examples 
#' s = simulate_MIMOSA2()
#' R = MIMOSA2(Ntot=s$Ntot, ns1 = s$ns1, nu1 = s$nu1, nu0 = s$nu0, ns0 = s$ns0)
#' Boxplot(obj = R,truth = s$truth)
Boxplot = function(obj,truth){
  with(obj,{
  ggplot(data.table(
    contrast = ps1_hat - pu1_hat - ps0_hat + pu0_hat,
    truth,
    inds = max.col(z)
  )) + geom_boxplot(outlier.colour = NA, aes(x = truth, y = contrast, fill =
                                               factor(inds))) + geom_jitter(
                                                 aes(x = truth, y = contrast, fill = factor(inds)),
                                                 position = position_jitterdodge(),
                                                 show.legend = FALSE
                                               ) + theme_bw() + scale_y_continuous("(ps1-pu1)-(ps0-pu0)") + scale_fill_discrete(
                                                 "Model Component",
                                                 labels = c(
                                                   "R1",
                                                   "NR",
                                                   "NSR1",
                                                   "NSR2",
                                                   "R2"
                                                 ),
                                                 limits = c(1:5)
                                               ) + scale_x_discrete("True class")
})
}

roc = function (p, truth) 
{
  s <- seq(0, 1, l = 1000)
  table <- t(sapply(s, function(th) {
    prop.table(table(test = factor(p <= th, levels = c("FALSE", 
                                                       "TRUE")), truth = truth), margin = 2)["TRUE", ]
  }))
  colnames(table) <- c("FPR", "TPR")
  table
}