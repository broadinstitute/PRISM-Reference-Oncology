
#' Fitting a dose response curve to the given dose and viability (FC) values.
#' This function fits 5 different dose-response functiosn to the given dose, viability pairs using dr4pl and drc packages and returns
#' the best one (lowest mse) among them.
#'
#' @param FC : Measured viability vector
#' @param dose : Dose vector corresponding to FC
#' @param UL_low : Lower limit for the upper asymptote
#' @param UL_up : Upper limit for the upper asympotote
#' @param slope_decreasing: Should the curve to be constrained to be decreasing or not.
#'
#' @return Returns a single row data-frame with following columns:
#'          fit_name : Name of the fitting method with the highest explained variance (lowest mse)
#'          lower_limit : Lower asmpytote for the selected fit
#'          upper_limit : Upper asymptote for the selected fit
#'          slope : The Hill slope for the selected fit
#'          inflection : inflection point of the selected fit (EC50)
#'          mse : mse of the selected fit
#'          mad : mad of the selected fit
#'          frac_var_explained : Explained variance for the best fit
#'          successful_fit: If any of the fitting methods has positive frac_var_explained
#'          auc_riemann : The average measured viability value across doses
#'          minimum_dose : Minimum measured dose
#'          maximum_dose : Maximum measured dose
#'          auc : auc of the log-dose vs viability curve normalized to the measured dose-range (takes values between 0 and 1)
#'          log2_ic50 : Log2 of IC50 for the fitted curve
#' @export
#'
#' @examples
get_best_fit <- function(FC, dose, UL_low=0.8, UL_up=1.01, slope_decreasing=TRUE) {
  require(dr4pl)
  require(drc)
  require(tidyverse)
  require(magrittr)
  
  
  # Fits a number of alternate models  to the DRC and chooses the best fit.
  
  # UL low is the lowerbound of UL we pass to the optimizer and UL_up is the upper bound of UL that we pass to the optimizer
  # fomat of output will be:-
  # results.df <- data.frame("fit_name"=character(),
  #                          "lower_limit"=double(),
  #                          "upper_limit"=double(),
  #                          "slope"=double(),
  #                          "inflection"=double(),
  #                          "mse"=double(), "mad" =double(),
  #                          "frac_var_explained"=double())
  
  dose = as.numeric(dose)
  FC = as.numeric(FC)
  slope_bound <- ifelse(slope_decreasing, 1e-5, Inf)  # bound the slopes by default unless passed another option
  riemann_auc <- mean(pmin(1,FC)) ## mean fold-change after rounding FC to 1.
  var_data = var(FC)
  
  minimum_dose = min(dose); maximum_dose = max(dose)
  
  results.df <- list(); ix = 1
  
  
  # FIT 1 ---
  drc_model <-  tryCatch(drc::drm(FC ~ dose, data= data.frame(FC = FC, dose = dose),
                                  fct=LL.4(names = c("slope", "Lower Limit", "Upper Limit", "ED50")),
                                  lowerl = c(-slope_bound,0.0, UL_low, -Inf),upperl = c(Inf,1.01,UL_up, Inf)),
                         error = function(e)
                         {return(list(convergence=FALSE, error=TRUE,fit=list(convergence=FALSE)))})
  # "slope" in drc package is -ve of slope in dr4pl package
  
  
  if (drc_model$fit$convergence & all(is.finite(unlist(drc_model$coefficients)))){
    mse_mad <- compute_mse_mad(FC, dose, as.numeric(drc_model$coefficients[[3]]), as.numeric(drc_model$coefficients[[2]]),
                               -as.numeric(drc_model$coefficients[[1]]), as.numeric(drc_model$coefficients[[4]]))
    # "slope" in drc package is -ve of slope in dr4pl package and so -ve sign needs to be put in here.
    
    results.df[[ix]] <- tibble( fit_name = "drc_drm_constrained",
                                lower_limit = as.numeric(drc_model$coefficients[[2]]),
                                upper_limit = as.numeric(drc_model$coefficients[[3]]),
                                slope = -as.numeric(drc_model$coefficients[[1]]),
                                inflection = as.numeric(drc_model$coefficients[[4]]),
                                mse = mse_mad$mse, mad = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
    ix = ix + 1
  }
  
  
  
  # FIT 2 ---
  drc_model <-  tryCatch(drc::drm(FC ~ dose, data= data.frame(FC = FC, dose = dose),
                                  fct=LL.4(names = c("slope", "Lower Limit", "Upper Limit", "ED50")),
                                  lowerl = c(-Inf,0.0, UL_low, -Inf),upperl = c(Inf,1.01,UL_up, Inf)),
                         error = function(e)
                         {return(list(convergence=FALSE, error=TRUE,fit=list(convergence=FALSE)))})
  # "slope" in drc package is -ve of slope in dr4pl package
  
  
  if (drc_model$fit$convergence & all(is.finite(unlist(drc_model$coefficients)))){
    if((!slope_decreasing) | (as.numeric(drc_model$coefficients[[1]]) > 0)){
      mse_mad <- compute_mse_mad(FC, dose, as.numeric(drc_model$coefficients[[3]]), as.numeric(drc_model$coefficients[[2]]),
                                 -as.numeric(drc_model$coefficients[[1]]), as.numeric(drc_model$coefficients[[4]]))
      # "slope" in drc package is -ve of slope in dr4pl package and so -ve sign needs to be put in here.
      
      results.df[[ix]] <- tibble( fit_name = "drc_drm_unconstrained",
                                  lower_limit = as.numeric(drc_model$coefficients[[2]]),
                                  upper_limit = as.numeric(drc_model$coefficients[[3]]),
                                  slope = -as.numeric(drc_model$coefficients[[1]]),
                                  inflection = as.numeric(drc_model$coefficients[[4]]),
                                  mse = mse_mad$mse, mad = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
      ix = ix + 1
    }
  }
  
  
  # FIT 3 ---
  dr4pl_initMan_optNM <- tryCatch(dr4pl(dose, FC,
                                        init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_2 = 8*min(dose),
                                                                       theta_3= -3, theta_4 = 0.01),
                                        lowerl = c(UL_low, -Inf, -Inf, 0),
                                        upperl = c(UL_up, Inf, slope_bound, 1.01),
                                        method.optim="Nelder-Mead"),
                                  error= function(e){return(list(convergence=FALSE, error=TRUE))}
  )
  
  if (!dr4pl_initMan_optNM$convergence){
    if (!is.null(dr4pl_initMan_optNM$dr4pl.robust)) {
      dr4pl_initMan_optNM <- dr4pl_initMan_optNM$dr4pl.robust
    }
  }
  
  if (dr4pl_initMan_optNM$convergence & all(is.finite(unlist(dr4pl_initMan_optNM$parameters)))){
    mse_mad <- compute_mse_mad(FC, dose, dr4pl_initMan_optNM$parameters[[1]], dr4pl_initMan_optNM$parameters[[4]],
                               dr4pl_initMan_optNM$parameters[[3]], dr4pl_initMan_optNM$parameters[[2]])
    
    results.df[[ix]] <- tibble( fit_name = "dr4pl_initMan_constrained_optNM",
                                lower_limit = as.numeric(dr4pl_initMan_optNM$parameters[[4]]),
                                upper_limit = as.numeric(dr4pl_initMan_optNM$parameters[[1]]),
                                slope = as.numeric(dr4pl_initMan_optNM$parameters[[3]]),
                                inflection = as.numeric(dr4pl_initMan_optNM$parameters[[2]]),
                                mse = mse_mad$mse, mad = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
    ix = ix + 1
  }
  
  # FIT 4 ---
  dr4pl_unconstrained <- tryCatch(dr4pl(dose, FC,
                                        init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.3),
                                        method.init = "logistic",
                                        lowerl = c(0.99, -Inf, -Inf, 0),
                                        upperl = c(1.01, Inf, Inf, 1.01)),
                                  error = function(e) {print(e); return(NA)})
  
  if (!all(is.na(dr4pl_unconstrained))) {
    if (!dr4pl_unconstrained$convergence) {
      dr4pl_unconstrained <- dr4pl_unconstrained$dr4pl.robust
    }
  }
  
  
  param <- tryCatch(dr4pl_unconstrained$parameters, error = function(e) return(NA))
  if (!all(is.na(param))){
    if(as.numeric(dr4pl_unconstrained$parameters[[3]])<slope_bound){ ### while slope bound is not passed to this last optimizer, we do not accept a solution not within the bound
      mse_mad <- compute_mse_mad(FC, dose, dr4pl_unconstrained$parameters[[1]], dr4pl_unconstrained$parameters[[4]],
                                 dr4pl_unconstrained$parameters[[3]], dr4pl_unconstrained$parameters[[2]])
      results.df[[ix]] <- tibble( fit_name = "dr4pl_initL_unconstrained",
                                  lower_limit = as.numeric(dr4pl_unconstrained$parameters[[4]]),
                                  upper_limit = as.numeric(dr4pl_unconstrained$parameters[[1]]),
                                  slope = as.numeric(dr4pl_unconstrained$parameters[[3]]),
                                  inflection = as.numeric(dr4pl_unconstrained$parameters[[2]]),
                                  mse = mse_mad$mse, mad = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
      ix = ix + 1
    }
  }
  
  # FIT 5 ---
  dr4pl_initL <- tryCatch(dr4pl(dose, FC,
                                init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.005),
                                method.init = "logistic",
                                lowerl = c(UL_low, -Inf, -Inf, 0),
                                upperl = c(UL_up, Inf, slope_bound, 1.01)),
                          error= function(e){return(list(convergence=FALSE, error=TRUE))}
  )
  
  if (!dr4pl_initL$convergence){
    if (!is.null(dr4pl_initL$dr4pl.robust)) {
      dr4pl_initL <- dr4pl_initL$dr4pl.robust
    }
  }
  
  if (dr4pl_initL$convergence & all(is.finite(unlist(dr4pl_initL$parameters)))){
    mse_mad <- compute_mse_mad(FC,dose, dr4pl_initL$parameters[[1]], dr4pl_initL$parameters[[4]],
                               dr4pl_initL$parameters[[3]], dr4pl_initL$parameters[[2]])
    
    results.df[[ix]] <- tibble( fit_name = "dr4pl_initL_constrained",
                                lower_limit = as.numeric(dr4pl_initL$parameters[[4]]),
                                upper_limit = as.numeric(dr4pl_initL$parameters[[1]]),
                                slope = as.numeric(dr4pl_initL$parameters[[3]]),
                                inflection = as.numeric(dr4pl_initL$parameters[[2]]),
                                mse = mse_mad$mse, mad = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
    ix = ix + 1
  }
  
  
  # Choose the best fit among the successful fits---
  results.df <- dplyr::bind_rows(results.df)
  
  if (nrow(results.df)>0){
    results.df <- results.df %>%
      dplyr::arrange(desc(frac_var_explained)) %>%
      head(1) %>%
      dplyr::mutate(successful_fit = TRUE,
                    auc_riemann = as.numeric(riemann_auc),
                    minimum_dose = minimum_dose, maximum_dose = maximum_dose) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(auc = compute_auc(lower_limit, upper_limit, inflection, slope, minimum_dose, maximum_dose),
                    log2_auc = log2(auc),
                    log2_ic50 = compute_log_ic50(lower_limit, upper_limit, inflection, slope, minimum_dose, maximum_dose))
    
  }else{
    results.df  <- data.frame(successful_fit=FALSE,
                              auc_riemann = as.numeric(riemann_auc),
                              minimum_dose = minimum_dose, maximum_dose = maximum_dose, auc = NA, log2_ic50 = NA)
  }
  
  return (results.df)
}




#' Fits a simple linear model by regressing y over each column of X.
#'
#' @param X : matrix of n by m, it can have NA's in it, columns should be named.
#' @param Y : matrix of n by p, rows are in the same order of with the rows of X.
#' @param v.th : Minimum variance for the columns of X, columns with smaller variances will be dropped.
#' @param n.min : Minimum number of finite pairs between columns of X and Y, column pairs not satisfying this condition will be dropped.
#' @param rank.max : For each column of Y, only the top rank.max elements will be returned.
#' @param q.val.max : Associations with a q_val > q.val.max will be dropped.
#' @param regression_coef : Should the regression coefficients be calculated
#' @param stability_score : Should the stability scores be calculated (valid only if regression_coef = TRUE) 
#' @param ns.min : If the stability scores are calculated, associations with score less than ns.min will be dropped
#' 
#'
#' @return A data frame with: x (corresponding column of X), y (corresponding column of Y), correlation_coeff (Pearson correlation), 
#'         regression_coeff (regression coefficient), p_val / q_val (homoskedastic p-value / q-value),
#'         n (number of non-na samples), ns (a stability score to represent how many points of y needs to be replaced with
#'         the mean value to flip the sign of regression coefficient), rank (rank of the significance for each column of Y), 
#
#' @export
#'
#' @examples
#'
linear_model <- function(X, Y, v.X.min = 0.0025, n.min = 25, rank.max = 250, q.val.max = 0.1, regression_coef = TRUE, stability_score = FALSE, ns.min = 3) {
  require(tidyverse)
  require(matrixStats)
  require(WGCNA)
  
  X <- X[, (matrixStats::colVars(X, na.rm = T) >= v.X.min) & (colSums(!is.na(X)) >= n.min), drop = FALSE]
  
  print("Calculating correlation coefficients and q-values...")
  res <- WGCNA::corAndPvalue(X,Y, use = "p")[c(1,2,5)] %>%
    reshape2::melt() %>% 
    tidyr::pivot_wider(names_from = L1, values_from = value) %>% 
    dplyr::filter(nObs >= n.min) %>% 
    dplyr::group_by(Var2) %>%  
    dplyr::mutate(q = p.adjust(p, method = "BH")) %>% 
    dplyr::arrange(q) %>% dplyr::mutate(rank = 1:n()) %>%
    dplyr::ungroup() %>% 
    dplyr::filter(q <= q.val.max) %>%
    dplyr::group_by(Var2, sign(cor)) %>% 
    dplyr::arrange(desc(abs(cor))) %>% dplyr::mutate(rank = pmin(rank, 1:n())) %>%
    dplyr::ungroup() %>% 
    dplyr::filter(rank <= rank.max) %>% 
    dplyr::rename(x = Var1, y = Var2, rho = cor, p.val = p, q.val = q, n = nObs) %>% 
    dplyr::distinct(x, y, rho, p.val, q.val, n, rank) %>%
    dplyr::rename(correlation_coef = rho, p_val = p.val, q_val = q.val)
  
  # Do we need the regression coefficients?
  if(regression_coef & (nrow(res) > 0)){
    print("Calculating regression coefficients...")
    X_ <- X[, res$x]; Y_ <- Y[, res$y]
    mask <- !is.finite(X_ * Y_)
    X_[mask] <- NA; Y_[mask] <- NA
    Y_ <- scale(Y_, center = TRUE, scale = FALSE)
    X_ <- scale(X_, center = TRUE, scale = FALSE)
    
    res$var.x <- colMeans(X_^2, na.rm = T)
    res$var.y <- colMeans(Y_^2, na.rm = T) 
    
    res <- res %>% 
      dplyr::mutate(regression_coef = correlation_coef * sqrt(var.y/var.x)) %>% 
      dplyr::select(-var.y)
    
    # Experimental robustness score?
    if(stability_score){
      print("Calculating stability scores...")
      S <- X_ * Y_
      n <- colSums(!is.na(X_))
      ns <- apply(S,2, function(s) sum(cumsum(sort(s/sum(s, na.rm = T), decreasing = T)) < 1, na.rm = T)) + 1
      res$ns <- pmin(ns, n - apply(X_, 2, function(x) max(table(x))))
      res <- dplyr::filter(res, ns >= ns.min)
    }
  }
  
  return(res)
}


#' Calculates and returns univariate analysis results by correlating each column of Y with features sets from depmap.org.
#'
#' @param Y : Matrix n x m, make sure rownames are depmap_id's and columns are named. There can be NA's.
#' @param file : Please point out to the downloaded depmap_datasets.h5 file.
#' @param features : You can give a subset of available feature sets, but if left as NULL, it computes for everything.
#' @param homoskedastic : Should homoskedastic or robust q-values be used for ranking and filtering. Default is homoskedastic.
#' @param n.X.min : Results with less thana given sample size are dropped, default is 100
#' @param n_stable.min : Results that can be flipped with less than n_stable.min data points are dropped, the default is 3
#' @param q_val.max : Results with q-values less than q_val.max are dropped, default is 0.2
#' @param parallel : Should multiple cores to be used to decrease the compute time. Default is FALSE.
#'
#' @return Returns a data-table with each row corresponds to a particular (feature, feature_set, y) triplet. See linear_model for the other columns.
#' @export
#'
#' @examples
#'
univariate_biomarker_table <- function(Y, file,
                                       features = NULL, n.X.min = 100,
                                       v.X.min = 0.0025,
                                       n_stable.min = 3, q_val_max = .1,
                                       regression_coef = TRUE,
                                       stability_score = TRUE, rank.max = 250){
  require(tidyverse)
  require(magrittr)
  require(rhdf5)
  require(WGCNA)
  
  
  if(!is.matrix(Y)){
    Y <- as.matrix(Y)
  }
  
  if(is.null(features)){
    features <-  substr(setdiff(h5ls(file)$group, "/"),2,100)
  }else{
    features <- intersect(features, substr(setdiff(h5ls(file)$group, "/"),2,100))
  }
  print(features)
  
  RESULTS <- list()
  for(feat in features){
    
    
    X <- read_dataset(file , feat)
    cl = intersect(rownames(X), rownames(Y))
    X <- X[cl, ]
    
    if((dim(X)[1] >= n.X.min) & (dim(X)[2] > 0)){
      
      print(paste0(feat, " - ", dim(X)[1] ,'x', dim(X)[2]))
      
      RESULTS[[feat]] <- linear_model(X = X, Y = Y[cl, , drop = FALSE],
                                      v.X.min = v.X.min, n.min = n.X.min, 
                                      rank.max = rank.max, q.val.max = q_val_max, 
                                      regression_coef = regression_coef, stability_score = stability_score, ns.min = n_stable.min) %>% 
        dplyr::rename(feature = x) %>%
        dplyr::mutate(feature_set = feat) 
    }
    print(feat)
  }
  
  RESULTS <- dplyr::bind_rows(RESULTS)
  return(RESULTS)
}

