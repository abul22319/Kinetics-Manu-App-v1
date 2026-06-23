packages = c("readxl", 
             "stats", 
             "ggplot2", 
             "lme4", 
             "dplyr", 
             "tidyverse", 
             "minpack.lm",
             "signal",
             "zoo")

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# Shared Helpers

R_Sqr = function(observed, predicted) {
  Residual_1  = sum((observed - predicted)^2)
  Residual_2  = observed - predicted
  Sum         = sum((observed - mean(observed))^2)
  R_squared   = 1 - (Residual_1 / Sum)
  return(list(R_squared = R_squared,
              Residuals = Residual_2))
}


filter_and_interpolate = function(Data, fc = 0.15) {
  fs      <- 1/3
  nyquist <- fs / 2
  Wn      <- fc / nyquist
  bf      <- butter(2, Wn, type = "low")
  
  Data$FBF     <- as.numeric(Data$FBF)
  Data$FBF     <- na.approx(Data$FBF, na.rm = FALSE)
  Data$FBF_filtered <- filtfilt(bf, Data$FBF)
  
  Data$FBF_rel <- as.numeric(Data$FBF_rel)
  Data$FBF_rel <- na.approx(Data$FBF_rel, na.rm = FALSE)
  Data$FBF_rel_filtered <- filtfilt(bf, Data$FBF_rel)
  
  Data$FVC     <- as.numeric(Data$FVC)
  Data$FVC     <- na.approx(Data$FVC, na.rm = FALSE)
  Data$FVC_filtered <- filtfilt(bf, Data$FVC)
  
  Data$FVC_rel <- as.numeric(Data$FVC_rel)
  Data$FVC_rel <- na.approx(Data$FVC_rel, na.rm = FALSE)
  Data$FVC_rel_filtered <- filtfilt(bf, Data$FVC_rel)
  
  return(Data)
}


make_plots = function(Data, col_filtered, col_fit, residuals) {
  
  exp_plot <- ggplot(Data, aes(x = Time)) +
    geom_line(aes(y = .data[[col_fit]]), linewidth = 0.8) +
    geom_point(aes(y = .data[[col_filtered]]), color = "blue", size = 2.5) +
    labs(x = "Time (s)", y = col_filtered) +
    theme_classic()
  
  ref_plot <- ggplot(Data, aes(x = Time)) +
    geom_point(aes(y = residuals), color = "darkorange", size = 2.5) +
    geom_hline(yintercept = mean(residuals), linewidth = 0.8) +
    labs(x = "Time (s)", y = "Residuals") +
    theme_classic()
  
  return(list(Exp.Model = exp_plot, RefLine.Model = ref_plot))
}

# Mono exp model

MonoExp = function(x, y, fc = 0.15) {
  
  Data <- read_excel(x, sheet = "Sheet1")
  Data <- filter_and_interpolate(Data, fc = fc)
  
  constrained <- (y == 1)
  
  fit_mono = function(col_name) {
    
    col_data <- Data[[col_name]]
    
    Start <- list(
      B   = 0.6 * (max(col_data) - min(col_data)) + min(col_data),
      tau = (0.2 * max(Data$Time)) / 3,
      TD1 = 0.1  * max(Data$Time)
    )
    
    formula_str <- as.formula(
      paste0(col_name,
             " ~ B * (1 - exp(-pmax(Time - TD1, 0) / tau))")
    )
    
    if (constrained) {
      Lo  <- c(B = 0, tau = 0, TD1 = 0)
      fit <- nlsLM(formula_str,
                   data    = Data,
                   start   = Start,
                   control = nls.lm.control(maxiter = 200),
                   lower   = Lo)
    } else {
      fit <- nlsLM(formula_str,
                   data    = Data,
                   start   = Start,
                   control = nls.lm.control(maxiter = 200))
    }
    
    cf <- coef(fit)
    
    params <- data.frame(
      B   = cf["B"],
      Tau = cf["tau"],
      TD  = cf["TD1"],
      MRT = cf["tau"] + cf["TD1"]
    )
    
    predicted <- predict(fit)
    rsqr_obj  <- R_Sqr(col_data, predicted)
    
    corr    <- cor.test(rsqr_obj$Residuals, col_data, method = "spearman")
    cor_df  <- data.frame(
      Rsq = rsqr_obj$R_squared,
      P   = corr$p.value,
      Rho = corr$estimate
    )
    
    col_fit          <- paste0(col_name, "_Fit")
    Data[[col_fit]]  <<- predicted
    
    plots <- make_plots(Data, col_name, col_fit, rsqr_obj$Residuals)
    
    return(list(
      Parameters    = params,
      Exp.Model     = plots$Exp.Model,
      RefLine.Model = plots$RefLine.Model,
      R.Sqr         = rsqr_obj$R_squared,
      Cor.Result    = cor_df
    ))
  }
  
  r_fbf     <- fit_mono("FBF_filtered")
  r_fbf_rel <- fit_mono("FBF_rel_filtered")
  r_fvc     <- fit_mono("FVC_filtered")
  r_fvc_rel <- fit_mono("FVC_rel_filtered")
  
  return(list(
    Parameters_FBF        = r_fbf$Parameters,
    Exp.Model_FBF         = r_fbf$Exp.Model,
    RefLine.Model_FBF     = r_fbf$RefLine.Model,
    R.Sqr_FBF             = r_fbf$R.Sqr,
    Cor.Result_FBF        = r_fbf$Cor.Result,
    
    Parameters_FBF_rel    = r_fbf_rel$Parameters,
    Exp.Model_FBF_rel     = r_fbf_rel$Exp.Model,
    RefLine.Model_FBF_rel = r_fbf_rel$RefLine.Model,
    R.Sqr_FBF_rel         = r_fbf_rel$R.Sqr,
    Cor.Result_FBF_rel    = r_fbf_rel$Cor.Result,
    
    Parameters_FVC        = r_fvc$Parameters,
    Exp.Model_FVC         = r_fvc$Exp.Model,
    RefLine.Model_FVC     = r_fvc$RefLine.Model,
    R.Sqr_FVC             = r_fvc$R.Sqr,
    Cor.Result_FVC        = r_fvc$Cor.Result,
    
    Parameters_FVC_rel    = r_fvc_rel$Parameters,
    Exp.Model_FVC_rel     = r_fvc_rel$Exp.Model,
    RefLine.Model_FVC_rel = r_fvc_rel$RefLine.Model,
    R.Sqr_FVC_rel         = r_fvc_rel$R.Sqr,
    Cor.Result_FVC_rel    = r_fvc_rel$Cor.Result
  ))
}

# Bi exp model
BiExp = function(x, y, fc = 0.15) {
  
  Data <- read_excel(x, sheet = "Sheet1")
  Data <- filter_and_interpolate(Data, fc = fc)
  
  constrained <- (y == 1)
  
  fit_bi = function(col_name) {
    
    col_data   <- Data[[col_name]]
    amp_range  <- max(col_data) - min(col_data)
    base       <- min(col_data)
    
    Start <- list(
      B1   = 0.4  * amp_range + base,
      tau1 = (0.1 * max(Data$Time)) / 3,
      TD1  = 0.05 * max(Data$Time),
      B2   = 0.2  * amp_range + base,
      tau2 = (0.4 * max(Data$Time)) / 3,
      TD2  = 0.1  * max(Data$Time)
    )
    
    formula_str <- as.formula(
      paste0(col_name,
             " ~ B1 * (1 - exp(-pmax(Time - TD1, 0) / tau1)) +",
             "   B2 * (1 - exp(-pmax(Time - TD2, 0) / tau2))")
    )
    
    if (constrained) {
      Lo  <- c(B1 = 0, tau1 = 0, TD1 = 0,
               B2 = 0, tau2 = 0, TD2 = 0)
      fit <- nlsLM(formula_str,
                   data    = Data,
                   start   = Start,
                   control = nls.lm.control(maxiter = 200),
                   lower   = Lo)
    } else {
      fit <- nlsLM(formula_str,
                   data    = Data,
                   start   = Start,
                   control = nls.lm.control(maxiter = 200))
    }
    
    cf <- coef(fit)
    
    params <- data.frame(
      B1   = cf["B1"],
      Tau1 = cf["tau1"],
      TD1  = cf["TD1"],
      MRT1 = cf["tau1"] + cf["TD1"],
      B2   = cf["B2"],
      Tau2 = cf["tau2"],
      TD2  = cf["TD2"],
      MRT2 = cf["tau2"] + cf["TD2"]
    )
    
    predicted <- predict(fit)
    rsqr_obj  <- R_Sqr(col_data, predicted)
    
    corr   <- cor.test(rsqr_obj$Residuals, col_data, method = "spearman")
    cor_df <- data.frame(
      Rsq = rsqr_obj$R_squared,
      P   = corr$p.value,
      Rho = corr$estimate
    )
    
    col_fit         <- paste0(col_name, "_Fit")
    Data[[col_fit]] <<- predicted
    
    plots <- make_plots(Data, col_name, col_fit, rsqr_obj$Residuals)
    
    return(list(
      Parameters    = params,
      Exp.Model     = plots$Exp.Model,
      RefLine.Model = plots$RefLine.Model,
      R.Sqr         = rsqr_obj$R_squared,
      Cor.Result    = cor_df
    ))
  }
  
  r_fbf     <- fit_bi("FBF_filtered")
  r_fbf_rel <- fit_bi("FBF_rel_filtered")
  r_fvc     <- fit_bi("FVC_filtered")
  r_fvc_rel <- fit_bi("FVC_rel_filtered")
  
  return(list(
    Parameters_FBF        = r_fbf$Parameters,
    Exp.Model_FBF         = r_fbf$Exp.Model,
    RefLine.Model_FBF     = r_fbf$RefLine.Model,
    R.Sqr_FBF             = r_fbf$R.Sqr,
    Cor.Result_FBF        = r_fbf$Cor.Result,
    
    Parameters_FBF_rel    = r_fbf_rel$Parameters,
    Exp.Model_FBF_rel     = r_fbf_rel$Exp.Model,
    RefLine.Model_FBF_rel = r_fbf_rel$RefLine.Model,
    R.Sqr_FBF_rel         = r_fbf_rel$R.Sqr,
    Cor.Result_FBF_rel    = r_fbf_rel$Cor.Result,
    
    Parameters_FVC        = r_fvc$Parameters,
    Exp.Model_FVC         = r_fvc$Exp.Model,
    RefLine.Model_FVC     = r_fvc$RefLine.Model,
    R.Sqr_FVC             = r_fvc$R.Sqr,
    Cor.Result_FVC        = r_fvc$Cor.Result,
    
    Parameters_FVC_rel    = r_fvc_rel$Parameters,
    Exp.Model_FVC_rel     = r_fvc_rel$Exp.Model,
    RefLine.Model_FVC_rel = r_fvc_rel$RefLine.Model,
    R.Sqr_FVC_rel         = r_fvc_rel$R.Sqr,
    Cor.Result_FVC_rel    = r_fvc_rel$Cor.Result
  ))
}
