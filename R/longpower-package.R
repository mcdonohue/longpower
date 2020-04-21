#' Sample size calculations for longitudinal data
#' 
#' The longpower package contains functions for computing power and sample size
#' for linear models of longitudinal data based on the formula due to Liu and
#' Liang (1997) and Diggle et al (1994). Either formula is expressed in terms
#' of marginal model or Generalized Estimating Equations (GEE) parameters. This
#' package contains functions which translate pilot mixed effect model
#' parameters (e.g.  random intercept and/or slope) into marginal model
#' parameters so that the formulas of Diggle et al or Liu and Liang formula can
#' be applied to produce sample size calculations for two sample longitudinal
#' designs assuming known variance. The package also handles the categorical
#' time Mixed Model of Repeated Measures (MMRM) using the formula of Lu, Luo,
#' and Chen (2008)
#' 
#' \tabular{ll}{ Package: \tab longpower\cr Type: \tab Package\cr Version: \tab
#' 1.0\cr Date: \tab 2013-05-22\cr License: \tab GPL (>= 2)\cr LazyLoad: \tab
#' yes\cr }
#' 
#' @name longpower-package
#' @docType package
#' @author Michael C. Donohue <mdonohue@@usc.edu> Anthony C. Gamst Steven D.
#' Edland
#' @seealso [`lmmpower``], [`power.mmrm`],
#' [`power.mmrm.ar1`], [`lmmpower`],
#' [`diggle.linear.power`], [`edland.linear.power`],
#' [`liu.liang.linear.power`]
#' @references Diggle PJ, Heagerty PJ, Liang K, Zeger SL. (2002) \emph{Analysis
#' of longitudinal data}. Second Edition.  Oxford Statistical Science Series.
#' 
#' Liu, G., & Liang, K. Y. (1997). Sample size calculations for studies with
#' correlated observations.  \emph{Biometrics}, 53(3), 937-47.
#' 
#' Lu, K., Luo, X., & Chen, P.-Y. (2008). Sample size estimation for repeated
#' measures analysis in randomized clinical trials with missing data.
#' \emph{International Journal of Biostatistics}, 4, (1)
#' @keywords package
NULL