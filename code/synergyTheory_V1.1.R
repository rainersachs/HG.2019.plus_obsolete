# Copyright:    (C) 2017-2019 Sachs Undergraduate Research Apprentice Program
#               This program and its accompanying materials are distributed 
#               under the terms of the GNU General Public License v3.
# Filename:     synergyTheory.R 
# Purpose:      Concerns radiogenic mouse Harderian gland tumorigenesis. 
#               Contains relevant synergy theory models, information coefficient
#               calculations and useful objects. It is part of the 
#               source code for the NASAmouseHG project.
# Contact:      Rainer K. Sachs 
# Website:      https://github.com/sachsURAP/NASAmouseHG
# Mod history:  18 Jun 2018
# Details:      See data_info.R for further licensing, attribution,
#               references, and abbreviation information.

source("data_info_V1.R") # Load in the data. 
# Remark: dose is in units of cGy; 
# LET usually in keV/micron; 
# prevalence Prev always < 1.
# (i.e. not in %, which would mean prevalence < 100 but is strongly deprecated).

library(deSolve) # Solving differential equations.
library(dplyr)
#========================= MISC. OBJECTS & VARIABLES ==========================#
# In next line phi controls how fast NTE build up from zero; not really needed 
# during calibration since phi * Dose >> 1 at every observed Dose !=0. phi is
# needed for later synergy calculations. 
# d_0 = 1 / phi = 5 x 10-4 cGy = 5 x 10^-6 Gy.

phi <- 2000 # even larger phi should give the same final results, 
            # but might cause extra problems with R. 
#Y_0 <- 0.046404
#Y_0 <- 0.025 # value for robustness check
#Y_0=0.041
print (Y_0)

#================================ DER MODELS ==================================#

#================ PHOTON MODEL =================#
# Linear model fit on beta_decay_data dataset. 
# We will never recalculate this unless new data comes in but here it is 
# just in case.

# beta_decay_lm <- lm(HG ~ dose + I(dose ^ 2), data = beta_decay_data)
# summary(beta_decay_lm, correlation = TRUE)

#=============== HZE/NTE MODEL =================#
# (HZE = high charge and energy; 
# select(filter(HZE_data, Beam == "O")
# NTE = non-targeted effects are included in addition to TE)
HZE_data <- select(filter(ion_data, Z > 3),1:length(ion_data[2,]))

# Uses 3 adjustable parameters. 
HZE_nte_model <- nls( # Calibrating parameters in a model that modifies the hazard function NTE models in 17Cuc. 
  Prev ~ Y_0 + (1 - exp ( - (aa1 * LET * dose * exp( - aa2 * LET) 
                               + (1 - exp( - phi * dose)) * kk1))), 
  data = HZE_data, 
  weights = NWeight,
  start = list(aa1 = .00009, aa2 = .001, kk1 = .06)) # Use extra argument trace = TRUE if you want to watch convergence. 

summary(HZE_nte_model, correlation = TRUE) # Parameter values & accuracy.
# If a paper uses dose in Gy some care is required in the preceding and following lines to rescale from cGy.
vcov(HZE_nte_model) # Variance-covariance matrix.
HZE_nte_model_coef <- coef(HZE_nte_model) # Calibrated central values of the 3 parameters.

# The DER, = 0 at dose 0.
calibrated_nte_hazard_func <- function(dose, LET, coef) { # Calibrated hazard function. 
 return(coef[1] * LET * dose * exp( - coef[2] * LET) 
        + (1 - exp( - phi * dose)) * coef[3])
} 

calibrated_HZE_nte_der <- function(dose, LET, coef = HZE_nte_model_coef) { # Calibrated HZE NTE DER.
  return(1 - exp( - calibrated_nte_hazard_func(dose, LET, coef)))
}

#================ HZE/TE MODEL =================#
# (TE = targeted effects only). This chunk runs and gives good results. 

HZE_te_model <- nls( # Calibrating parameters in a TE only model.
  Prev ~ Y_0 + (1 - exp ( - (aate1 * LET * dose * exp( - aate2 * LET)))),
  data = HZE_data,
  weights = NWeight,
  start = list(aate1 = .00009, aate2 = .01))

summary(HZE_te_model, correlation = TRUE) # Parameter values & accuracy.
vcov(HZE_te_model) # Variance-covariance matrix.
HZE_te_model_coef <- coef(HZE_te_model) # Calibrated central values of the 2 parameters. 

# The DER, = 0 at dose 0.
calibrated_te_hazard_func <- function(dose, LET, coef) { # Calibrated hazard function.
  return(coef[1] * LET * dose * exp( - coef[2] * LET))
}

calibrated_HZE_te_der <- function(dose, LET, coef = HZE_te_model_coef) {
  return(1 - exp( - calibrated_te_hazard_func(dose, LET, coef))) # Calibrated HZE TE one-ion DER.
}

#==== LIGHT ION, LOW Z (<= 3), LOW LET MODEL ===#
low_LET_data = ion_data[c(1:12,48), ] # Swift protons and alpha particles. RKS to RKS: need to add extra proton point
low_LET_model <- nls(
  Prev ~ Y_0 + 1 - exp( - alpha_low * dose), # alpha is used throughout radiobiology for dose coefficients.
  data = low_LET_data,
  weights = NWeight,
  start = list(alpha_low = .005))

summary(low_LET_model, correlation = TRUE)
low_LET_model_coef <- coef(low_LET_model) # Calibrated central values of the parameter.

# Calibrated Low LET model. Use L = 0, but maybe later will use small L > 0.
calibrated_low_LET_der <- function(dose, LET, alph_low = low_LET_model_coef[1]) {  
  return(1 - exp( - alph_low * dose))
}  

# Slope dE/dd of the low LET, low Z model.
low_LET_slope <- function(dose, LET) { 
  low_LET_model_coef * exp( - low_LET_model_coef * dose)  
}


#=========================== INFORMATION CRITERION ============================#
info_crit_table <- cbind(AIC(HZE_te_model, HZE_nte_model), 
                         BIC(HZE_te_model, HZE_nte_model))
print(info_crit_table)

##=================== Cross validation ====================##
## ====== ======##
# Seperate Data into 8 blocks, i.e. test/training sets:
data_len <- 1:length(HZE_data)
O_350 <- select(filter(HZE_data, LET == 20), data_len) #RKS added 3 oxygen points 5/17/2019)
Ne_670 <- select(filter(HZE_data, LET == 25), data_len) # Beam data not labeled
Si_260 <- select(filter(HZE_data, LET == 70), data_len)
Ti_1000 <- select(filter(HZE_data, LET == 100), data_len)
Fe_600 <- select(filter(HZE_data, LET == 193), data_len)
Fe_350 <- select(filter(HZE_data, LET == 253), data_len)
Nb_600 <- select(filter(HZE_data, LET == 464), data_len)
La_593 <- select(filter(HZE_data, LET == 953), data_len)

set_list <- list(O_350, Ne_670, Si_260, Ti_1000, Fe_600, Fe_350, Nb_600, La_593)
actual_prev <- HZE_data$Prev

# # Cross Validation for NTE Model:
theoretical <- vector()

for (i in 1:length(set_list)) {
  test <- set_list[[i]]
  excluded_list <- set_list[-i]
  train <- excluded_list[[1]]
  
  for (j in 2:length(excluded_list)) {
    train <- rbind(train, excluded_list[[j]])
  }
  
  HZE_nte_model <- nls(Prev ~ Y_0 + (1 - exp ( - (aa1 * LET * dose * exp( - aa2 * LET) + (1 - exp( - phi * dose)) * kk1))),
                       data = train,
                       weights = NWeight,
                       start = list(aa1 = .00009, aa2 = .001, kk1 = .06))
  predic <- predict(HZE_nte_model, test)
  theoretical <- c(theoretical, predic)
}

errors <- (theoretical - actual_prev)^2
NTE_cv <- weighted.mean(errors, HZE_data$NWeight)

#======= Cross Validation for TE Model  ========#
theoretical <- vector()
for (i in 1:8) {
  test <- set_list[[i]]
  excluded_list <- set_list[-i]
  train <- excluded_list[[1]]
  
  for (j in 2:length(excluded_list)) {
    train <- rbind(train, excluded_list[[j]])
  }
  
  HZE_te_model <- nls(Prev ~ Y_0 + (1 - exp ( - (aate1 * LET * dose * exp( - aate2 * LET)))),
                      data = train,
                      weights = NWeight,
                      start = list(aate1 = .00009, aate2 = .01))
  predic <- predict(HZE_te_model, test)
  theoretical <- c(theoretical, predic)
}

errors <- (theoretical - actual_prev)^2
TE_cv <- weighted.mean(errors, HZE_data$NWeight)

CV_table <- cbind(NTE_cv, TE_cv)
CV_table
# Dataframe subsetting ---------------------------------------------------------

for (i in 1:length(HZE_data)) { # Columns by index
  print(HZE_data[, i])
}

for (i in colnames(HZE_data)) { # Columns by names
  print(HZE_data[, i])
}

for (i in 1:length(HZE_data[1])) { # Rows by index
  print(HZE_data[i, ])
}

for (i in rownames(HZE_data)) { # Rows by name
  print(HZE_data[i, ])
}

# More sophisticated methods here: 
# https://www.r-bloggers.com/5-ways-to-subset-a-data-frame-in-r/


## == Rprojects/Auxiliary_files gives toy calculation of completion == ##
## ======= also see, e.g., R_projects/Old_projects/Chang_2019 ======== ##
## === But there it is inextricably tangled up with CRAN examples ==== ##


## ===== BASELINE NO-SYNERGY/ANTAGONISM MIXTURE DER FUNCTIONS ======== ##

#======== SIMPLE EFFECT ADDITIVITY (SEA) =======#

#' @description Applies Simple Effect Additivity to get a baseline mixture DER.
#' 
#' @param dose Numeric vector corresponding to the sum dose in cGy.
#' @param LET Numeric vector of all LET values, must be length n.
#' @param ratios Numeric vector of all dose ratios, must be length n.
#' @param lowLET Boolean of whether a low LET DER should be included in the mixture DER.
#' @param n Number of DERs, optional argument used to check parameter validity.
#' 
#' @details Corresponding elements of ratios, LET should be associated with the
#'          same DER.
#'          
#' @return Numeric vector representing the estimated Harderian Gland 
#'         prevalence from a SEA mixture DER constructed from the given DER 
#'         parameters. 
#'         
#' @examples
#  calculate_SEA(.01 * 0:40, c(70, 195), c(1/2, 1/2), n = 2)
#  calculate_SEA(.01 * 0:70, c(0.4, 195), c(4/7, 3/7))

calculate_SEA <- function(dose, LET, ratios, lowLET = FALSE, n = NULL) {
  if (!is.null(n) && (n != length(ratios) | n != length(LET))) {
    stop("Length of arguments do not match.") 
  } else if (sum(ratios) != 1) {
    stop("Sum of ratios do not add up to one.")
  } #  End error handling
  
  total <- 0
  i <- 1
  
  if (lowLET == TRUE) { 
    # First elements of ratios and LET should be the low-LET DER
    total <- total + calibrated_low_LET_der(dose * ratios[i], LET[i])
    i <- i + 1
  } 
  
  while (i < length(ratios) + 1) { # Iterate over HZE ions in the mixture.
    total <- total + calibrated_HZE_nte_der(dose * ratios[i], LET[i])
    i <- i + 1
  }
  
  return(total)
}


#====== INCREMENTAL EFFECT ADDITIVITY (IEA) ====#

#' @description Applies IEA to get a baseline no-synergy/antagonism mixture DER.
#' 
#' @param dose Numeric vector corresponding to the sum dose in cGy.
#' @param LET Numeric vector of all LET values, must be length n.
#' @param ratios Numeric vector of all dose ratios, must be length n.
#' @param model String value corresponding to the model to be used, either "NTE" or "TE". 
#' @param coef Named list of numeric vectors containing coefficients for one-ion DERs.
#' @param ders Named list of functions containing relevant one-ion DER models.
#' @param calculate_dI Named vector of functions to calculate dI depending on 
#'                     the selected model.
#' @param phi Numeric value, used in NTE models.
#' 
#' @details Corresponding elements of ratios, LET should be associated with the
#'          same DER.
#'          
#' @return Numeric vector representing the estimated Harderian Gland 
#'         prevalence from an IEA mixture DER constructed from the given 
#'         one-ion DERs parameters. 
#'         
#' @examples
#  calculate_id(.01 * 0:40, c(70, 195), c(1/2, 1/2))
#  calculate_id(.01 * 0:70, c(.4, 195), c(4/7, 3/7), model = "TE")

calculate_id <- function(dose, LET, ratios, model = "NTE",
                         coef = list(NTE = HZE_nte_model_coef, 
                                     TE = HZE_te_model_coef, 
                                     lowLET = low_LET_model_coef),
                         ders = list(NTE = calibrated_HZE_nte_der, 
                                     TE = calibrated_HZE_te_der, 
                                     lowLET = calibrated_low_LET_der),
                         calculate_dI = c(NTE = .calculate_dI_nte, 
                                          TE = .calculate_dI_te),
                         phi = 2000) {
  dE <- function(yini, state, pars) { # Constructing an ODE from the DERS.
    with(as.list(c(state, pars)), {

      # Screen out low LET values.
      lowLET_total <- lowLET_ratio <- 0
      remove <- c()
      for (i in 1:length(LET)) {
        if (LET[i] <= 3) { # Ion is low-LET.
          lowLET_total <- lowLET_total + LET[i]
          lowLET_ratio <- lowLET_ratio + ratios[i]
          remove <- unique(c(remove, LET[i]))
          ratios[i] <- 0
        }
      }
      LET <- c(setdiff(LET, remove))
      ratios <- ratios[! ratios == 0]
      
      # Begin calculating dI values.
      aa <- u <- dI <- vector(length = length(LET))
      if (length(LET) > 0) {
        for (i in 1:length(LET)) { 
          aa[i] <- pars[1] * LET[i] * exp( - pars[2] * LET[i])
          u[i] <- uniroot(function(dose) HZE_der(dose, LET[i], pars) - I,  
                          interval = c(0, 20000), 
                          extendInt = "yes", 
                          tol = 10 ^ - 10)$root 
          dI[i] <- ratios[i] * calc_dI(aa[i], u[i], pars[3]) 
        }
      }
      if (lowLET_ratio > 0) { 
        # If low-LET DER is present then include it at the end of the dI vector.
        u[length(LET) + 1] <- uniroot(function(dose) 
                                      low_der(dose, LET = lowLET_total, 
                                      alph_low = coef[["lowLET"]]) - I, 
                                      interval = c(0, 20000), 
                                      extendInt = "yes", 
                                      tol = 10 ^ - 10)$root
        dI[length(LET) + 1] <- lowLET_ratio * low_LET_slope(u[length(LET) + 1], 
                                                            LET = lowLET_total)
      }
      return(list(sum(dI)))
    }) 
  }
  p <- list(pars = coef[[model]], 
            HZE_der = ders[[model]], 
            low_der = ders[["lowLET"]], 
            calc_dI = calculate_dI[[model]])
  return(ode(c(I = 0), times = dose, dE, parms = p, method = "radau"))
}

#============= dI HIDDEN FUNCTIONS =============#
.calculate_dI_nte <- function(aa, u, kk1) {
  return((aa + exp( - phi * u) * kk1 * phi) * 
           exp( - (aa * u + (1 -exp( - phi * u)) * kk1)))
}

.calculate_dI_te <- function(aa, u, pars = NULL) {
  return(aa * exp(- aa * u))
}

#============================ DEVELOPER FUNCTIONS =============================#
test_runtime <- function(f, ...) { # Naive runtime check
  start_time <- Sys.time()
  f(...)
  Sys.time() - start_time
}
