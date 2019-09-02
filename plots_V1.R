# Copyright:    (C) 2017-2019 Sachs Undergraduate Research Apprentice Program
#               This program and its accompanying materials are distributed 
#               under the terms of the GNU General Public License v3.
# Filename:     plots.R 
# Purpose:      Concerns radiogenic mouse Harderian gland tumorigenesis. 
#               Contains code to generate many of the LSSR paper figures. It's
#               part of the source code for the Chang 2019 HG project.
# Contact:      Rainer K. Sachs 
# Website:      https://github.com/rainersachs/mouseHG_Chang_2019plus
# Mod history:  04 Apr 2019 # Edward: some comments in July 2019
# Details:      See data_info.R for further licensing, attribution, 
#               references, and abbreviation information.

source("monteCarlo_V1.R") # Load Monte Carlo.
library(Hmisc) # Error bars.
# library(dplyr)# already called in synergy

#= shift needed due to the synergy-theory-oriented way DERs are calibrated =#
dat_down = ion_data 
# dat_down stands for ion data shifted downward (by subtracting background) 
vvector = ion_data[ , "Prev"] # vvector is used only in this chunk
vvector=vvector - Y_0
dat_down[ , "Prev"]= vvector
# Edward: above chunk added 6/2/2019 

# ==========================================================#
#====== Fig. 1_final. Convex, Concave, Standard ============#
#=== Figs. 1&2_final are schematic. need not be redrawn ====#
# ==========================================================#
# d2 <- 0.01 * 0:200
# a <- 2; b <- .6; c <- 4
# E1 <- a * d2 + b * d2 ^ 2  # Convex
# E2 <- a * d2  # Linear no-threshold (LNT) same initial slope
# E3 <- c * (1 - (exp(- a * d2 / c))) # Concave, same initial slope
# plot(d2, E1, type = 'l', lwd = 3, bty = 'l', ann = FALSE)
# lines(d2, E2, lwd = 3)
# lines(d2, E3, lwd = 3)
# 
# a <- 0.45; b <- 1 / 8; c <- 0.8
# E1 <- b * d2 + 0.35 * b * d2 ^ 2  # Convex
# E2 <- 0.7 * a * d2  # Linear no-threshold (LNT) 
# E3 <- c * (1 - (exp(- 2 * a * d2 / c))) # Concave
# plot(d2, E3, type = 'l', lwd = 3, bty = 'l', ann = FALSE)
# lines(d2, E1, lwd = 3)
# lines(d2, E2, lwd = 3)
#============================================================#
#= Fig.3_final. Low LET Data SHIFTED, Error Bars, and DER. ==#
#============================================================#
ddose <- 0:701 
plot(c(0, 701), c(-.02, .75), pch = 19, col = 'white', ann = FALSE, bty = 'u') #RKS to RKS: force 7 ticks on y axis, see par arguments to plot()
lines(ddose, 1 - exp(- coef(summary(low_LET_model, correlation = TRUE))[1] * ddose), lwd = 2)# next is alpha particle data

errbar(dat_down[dat_down$Z == 2, ][, "dose"], 
       dat_down[dat_down$Z == 2, ][, "Prev"],
       yplus  = dat_down[dat_down$Z == 2, ][, "Prev"] + dat_down[
         dat_down$Z == 2, ][, "SD"], 
       yminus = dat_down[dat_down$Z == 2, ][, "Prev"] - dat_down[
         dat_down$Z == 2, ][, "SD"], pch = 19, cap = 0.02, add = TRUE,
       col = 'red', errbar.col = 'red', lwd = 1) # Proton data
legend(x = "bottomright", legend = "SD", cex=0.6)
print(Y_0) # amount data points shifted down compared to DER curves etc.
# Alpha particle data

errbar(dat_down[dat_down$Z == 1, ][, "dose"], 
       dat_down[dat_down$Z == 1, ][, "Prev"],
       yplus  = dat_down[dat_down$Z == 1, ][, "Prev"] + dat_down[
         dat_down$Z == 1, ][, "SD"], 
       yminus = dat_down[dat_down$Z == 1, ][, "Prev"] - dat_down[
         dat_down$Z == 1, ][, "SD"], pch = 19, cap = 0.02, add = TRUE,
       col = 'black', errbar.col = 'black', lwd = 1) # Proton data
legend(x = "bottomright", legend = "SD", cex=0.6)
print(Y_0) # amount data points shifted down compared to DER curves etc.

##==================================================#
#== Fig. for ICRR August 2019 Fe DERs. ==#
##==================================================#
d1_Fe <- c(0.001 * 0:9, 0.01 * 1:9,.1*1:9, 1:20)
fe_six_1 = calibrated_HZE_te_der(dose = d1_Fe, L = 193) #TE-only Fe600 1-ion DER
plot(c(0, 20.1), c(0, .201), col = "white", bty = 'L', ann = FALSE) # Set plot area
lines(d1_Fe, fe_six_1, col = 'blue') # Fe TE-only 1-ion DER, no background
fe_six1_nte = calibrated_HZE_nte_der(dose = d1_Fe, L = 193) # same for NTE-also
lines(d1_Fe, fe_six1_nte, col = 'black', lwd =2)

##==================================================#
#== Fig.4A_final Fe DERs. points, error bars, SHIFTED==#
##==================================================#
d1_Fe <- c(0.001 * 0:9, 0.01 * 1:9,.1*1:9, 1:80)
fe_six = calibrated_HZE_te_der(dose = d1_Fe, L = 193) #TE-only Fe600 1-ion DER
plot(c(0, 80.1), c(0, .65), col = "white", bty = 'L', ann = FALSE) # Set plot area
lines(d1_Fe, fe_six, col = 'black') # Fe TE-only 1-ion DER, no background
fe_six_nte = calibrated_HZE_nte_der(dose = d1_Fe, L = 193) # same for NTE-also
lines(d1_Fe, fe_six_nte, col = 'red', lwd =2)

errbar(dat_down[dat_down$LET==193,][,"dose"], dat_down[dat_down$LET==193,][,"Prev"],
       yplus  = dat_down[dat_down$LET == 193, ][, "Prev"] + dat_down[dat_down$LET == 193, ][, "SD"], 
       yminus = dat_down[dat_down$LET == 193, ][, "Prev"] - dat_down[dat_down$LET == 193, ][, "SD"],
       pch = 19, cap = 0.04, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1)
legend(x = "topleft", legend = "Fe 600, not 95%CI, SD, Y_0=", cex=0.49)
print(Y_0)

## =============================================================#
##= Fig4_B final: Fe 1-ion NTE-also der "mixture" + both ribbons =# 
##==============================================================#
# Declare ratios and LET values for plot; Fe with LET 193 
dd0 <- c(0.01 * 0:9, 0.1 * 1:9, 1:81)

# We use the plot that neglects adjustable parameter correlations
uncorr_fig_0 <- simulate_monte_carlo(n = 500, dd0, 193, 1, vcov = FALSE)

ci_data <- data.frame(dose = dd0,
                      # Monte Carlo values
                      uncorrBottom = uncorr_fig_0$monte_carlo[1, ],
                      uncorrTop = uncorr_fig_0$monte_carlo[2, ], 
                      
                      # one-ion DER for comparison
                      fe_six = calibrated_HZE_nte_der(dose = dd0, L = 193),
                      # IEA baseline mixture DER I(d), denoted by id below
                      i = calculate_id(dd0, 193, 1, model = "NTE")[, 2])

plot(c(0, 81), c(0, .62), col = "white", bty = 'L', ann = FALSE) # Set plot area

polygon(x = c(dd0, rev(dd0)), 
        y = c(ci_data[, "uncorrTop"], rev(ci_data[, "uncorrBottom"])), 
        col = "aquamarine2", lwd = .4, border = "aquamarine2") # CI ribbon

# We use the plot that takes adjustable parameter correlations into account
corr_fig_0 <- simulate_monte_carlo(n = 500, dd0, 193, 1)
ci_data <- data.frame(dose = dd0,
                      # Monte Carlo values
                      corrBottom = corr_fig_0$monte_carlo[1, ],
                      corrTop = corr_fig_0$monte_carlo[2, ], #
                      
                      # one-ion DERs for comparison
                      fe_six = calibrated_HZE_nte_der(dose = dd0, L = 193),
                      # p = calibrated_low_LET_der(dose = dd0, L = .4),
                      
                      # IEA baseline mixture DER I(d), denoted by id below
                      i = calculate_id(dd0, 193, 1, model = "NTE")[, 2])

polygon(x = c(dd0, rev(dd0)), 
        y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        col = "yellow", lwd = .4, border = "orange") # CI ribbon
lines(dd0, calibrated_HZE_nte_der(dose = dd0, L = 193), col = 'red', lwd = 2)
# ========================================================== #
#= Fig. 5A_final p-Fe mix, NTE, points+error, narrow ribbon =# 
# ========================================================== #
# Declare ratios and LET values for plot
ratios <- c(3/7, 4/7) # for Fe-p 
LET_vals <- c(193, .4)
d5_new <- c(0.01 * 0:9, 0.1 * 1:9, 1:70)
d5_Fe = c(0.01 * 0:9, 0.1 * 1:9, 1:30)
d5_p <- c(0.01 * 0:9, 0.1 * 1:9, 1:40)
# We use the plot that takes adjustable parameter correlations into account
corr_fig_5A <- simulate_monte_carlo(n = 500, d5_new, LET_vals, ratios, model = "NTE")
# The first argument, n, is the number of Monte Carlo repeats. Increase for
# greater accuracy. Decrease to speed up the program.
ci_data <- data.frame(dose = d5_new,
                      # Monte Carlo values
                      corrBottom = corr_fig_5A$monte_carlo[1, ],
                      corrTop = corr_fig_5A$monte_carlo[2, ], #
                      
                      # one-ion DERs 
                      fe_six = calibrated_HZE_nte_der(dose = d5_new, L = 193),
                      p = calibrated_low_LET_der(dose = d5_new, L = .4),
                      
                      # IEA baseline mixture DER I(d), denoted by id below
                      i = calculate_id(d5_new, LET_vals, ratios, model = "NTE")[, 2])

plot(c(0, 70.1), c(0, 0.601), col = "white", bty = 'l', ann = FALSE) # Set plot area
polygon(x = c(d5_new, rev(d5_new)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", lwd = .4, border = "orange") # Narrow CI ribbon
fe_six_one = calibrated_HZE_nte_der(dose = d5_Fe, L = 193)
p_one = calibrated_low_LET_der(dose = d5_p, L = .4)
lines(d5_p, p_one, col = 'brown', lwd = 2) # p DER
lines(d5_Fe, fe_six_one, col = 'blue') # Fe NTE-also 1-ion DER
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 1) # I(d)

Prev=as.numeric(mix_data$Prev)-Y_0; SD=as.numeric(mix_data$SD)
errbar(70,Prev[3], yplus  = Prev[3] +  SD[3],
       yminus = Prev[3] -SD[3],
       pch = 19, cap = 0.05, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1)
#legend(x = "topleft", legend = "not 95%CI, SD", cex=0.3)

# ========================================================== #
#= Fig. 5B_final p-Fe mix, TE, point & error, narrow ribbon =# 
# ========================================================== #
# We use the plot that takes adjustable parameter correlations into account
corr_fig_5B <- simulate_monte_carlo(n = 500, d5_new, LET_vals, ratios, model = "TE")
# The first argument, n, is the number of Monte Carlo repeats. Increase for
# greater accuracy. Decrease to speed up the program.
ci_data <- data.frame(dose = d5_new,
                      # Monte Carlo values
                      corrBottom = corr_fig_5B$monte_carlo[1, ],
                      corrTop = corr_fig_5B$monte_carlo[2, ], #
                      
                      # one-ion DERs 
                      fe_six = calibrated_HZE_te_der(dose = d5_new, L = 193),
                      p = calibrated_low_LET_der(dose = d5_new, L = .4),
                      
                      # IEA baseline mixture DER I(d), denoted by id below
                      i = calculate_id(d5_new, LET_vals, ratios, model = "TE")[, 2])

plot(c(0, 70.1), c(0, 0.601), col = "white", bty = 'l', ann = FALSE) # Set plot area
polygon(x = c(d5_new, rev(d5_new)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", lwd = .4, border = "orange") # Narrow CI ribbon
fe_six_te_one = calibrated_HZE_te_der(dose = d5_Fe, L = 193)
p_one = calibrated_low_LET_der(dose = d5_p, L = .4)
lines(d5_p, p_one, col = 'brown', lwd = 2) # p DER
lines(d5_Fe, fe_six_te_one, col = 'blue') # Fe TE-only 1-ion DER
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 1) # I(d)


Prev=as.numeric(mix_data$Prev)-Y_0; SD=as.numeric(mix_data$SD)
errbar(70,Prev[3], yplus  = Prev[3] +  SD[3],
       yminus = Prev[3] -SD[3],
       pch = 19, cap = 0.05, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1)


# ================================================================ #
#=== Figs. 6&7_A&B_final. 2 other mixtures. NTE-also, TE-only. ====#
#== Export as images width 400, height 350, font arial. then jpg ==#
#= Fig.6A_final. Fe-Si 50-50 total 40 cGy, NTE-also. like Fig. 5A =#
#==================================================================#
# Declare ratios and LET values for plot
ratios <- c(1/2, 1/2)
LET_vals <- c(193, 70)
d6 <- c(0.01 * 0:9, 0.1 * 1:9, 1:40)
d6_ion <- c(0.01 * 0:9, 0.1 * 1:9, 1:20)
# We use the plot that takes adjustable parameter correlations into account
corr_fig_6A <- simulate_monte_carlo(n = 500, d6, LET_vals, ratios, model = "NTE")
# The first argument, n, is the number of Monte Carlo repeats. Increase for
# greater accuracy. Decrease to speed up the program.
ci_data <- data.frame(dose = d6,
                      # Monte Carlo values
                      corrBottom = corr_fig_6A$monte_carlo[1, ],
                      corrTop = corr_fig_6A$monte_carlo[2, ], #
                      
                      # one-ion DERs for comparison
                      fe_six = calibrated_HZE_nte_der(dose = d6, L = 193),
                      si = calibrated_HZE_nte_der(dose = d6, L = 70),
                      # does d6_ion above rather than d_6 work?
                      # IEA baseline mixture DER I(d), denoted by id below
                      i = calculate_id(d6, LET_vals, ratios, model = "NTE")[, 2])

#We make the ribbon plot for correlated parameters
plot(c(0, 41), c(0, .36), col = "white", bty = 'L', ann = FALSE) # Set plot area
polygon(x = c(d6, rev(d6)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", lwd = .4, border = "orange") # Narrow CI ribbon

fe_six_one = calibrated_HZE_nte_der(dose = d6_ion, L = 193)
si_one = calibrated_HZE_nte_der(dose = d6_ion, L = 70)
lines(d6_ion, si_one, col = 'brown', lwd=2)
lines(d6_ion, fe_six_one, col = 'blue', lwd=2) # Fe TE-only 1-ion DER
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 2) # IEA NSNA I(d)

Prev=as.numeric(mix_data$Prev)-Y_0; SD=as.numeric(mix_data$SD)
errbar(40,Prev[5], yplus  = Prev[5] +  SD[5],
       yminus = Prev[5] -SD[5],
       pch = 19, cap = 0.05, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1)

#=========== 6B_final. repeat Fig. 6A_final but for TE-only DERs ============#
# We use the plot that takes adjustable parameter correlations into account
corr_fig_6B <- simulate_monte_carlo(n = 500, d6, LET_vals, ratios, model = "TE")
# The first argument, n, is the number of Monte Carlo repeats. Increase for
# greater accuracy. Decrease to speed up the program.
ci_data <- data.frame(dose = d6,
                      # Monte Carlo values
                      corrBottom = corr_fig_6B$monte_carlo[1, ],
                      corrTop = corr_fig_6B$monte_carlo[2, ], #
                      
                      # one-ion DERs for comparison
                      fe_six = calibrated_HZE_te_der(dose = d6, L = 193),
                      si = calibrated_HZE_te_der(dose = d6, L = 70),
                      # does d6_ion above rather than d_6 work?
                      # IEA baseline mixture DER I(d), denoted by id below
                      i = calculate_id(d6, LET_vals, ratios, model = "TE")[, 2])

#We make the ribbon plot for correlated parameters
plot(c(0, 41), c(0, .36), col = "white", bty = 'L', ann = FALSE) # Set plot area
polygon(x = c(d6, rev(d6)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", lwd = .4, border = "orange") # Narrow CI ribbon

fe_six_one = calibrated_HZE_te_der(dose = d6_ion, L = 193)
si_one = calibrated_HZE_te_der(dose = d6_ion, L = 70)
lines(d6_ion, si_one, col = 'brown', lwd=2)
lines(d6_ion, fe_six_one, col = 'blue', lwd=2) # Fe TE-only 1-ion DER
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 2) # IEA NSNA I(d)

Prev=as.numeric(mix_data$Prev)-Y_0; SD=as.numeric(mix_data$SD)
errbar(40,Prev[6], yplus  = Prev[6] +  SD[6],
       yminus = Prev[6] -SD[6],
       pch = 19, cap = 0.05, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1)

# ========================================================== #
#= Fig. 7A_final p-Si mix, NTE, points+error, narrow ribbon =# 
# ========================================================== #
# Declare ratios and LET values for plot
ratios <- c(.6, .4) # for Fe-p 
LET_vals <- c(.4, 70)
d7_new <- c(0.01 * 0:9, 0.1 * 1:9, 1:100)
d7_Si = c(0.01 * 0:9, 0.1 * 1:9, 1:40)
d7_p <- c(0.01 * 0:9, 0.1 * 1:9, 1:60)
# We use the plot that takes adjustable parameter correlations into account
corr_fig_7A <- simulate_monte_carlo(n = 500, d7_new, LET_vals, ratios, model = "NTE")
# The first argument, n, is the number of Monte Carlo repeats. Increase for
# greater accuracy. Decrease to speed up the program.
ci_data <- data.frame(dose = d7_new,
                      # Monte Carlo values
                      corrBottom = corr_fig_7A$monte_carlo[1, ],
                      corrTop = corr_fig_7A$monte_carlo[2, ], #
                      
                      # one-ion DERs 
                      Si = calibrated_HZE_nte_der(dose = d7_new, L = 70),
                      p = calibrated_low_LET_der(dose = d7_new, L = .4),
                      
                      # IEA baseline mixture DER I(d), denoted by id below
                      i = calculate_id(d7_new, LET_vals, ratios, model = "NTE")[, 2])

plot(c(0, 101), c(0,.36), col = "white", bty = 'L', ann = FALSE) # Set plot area
polygon(x = c(d7_new, rev(d7_new)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", lwd = .4, border = "orange") # Narrow CI ribbon
si_one = calibrated_HZE_nte_der(dose = d7_Si, L = 70)
p_one = calibrated_low_LET_der(dose = d7_p, L = .4)
lines(d7_p, p_one, col = 'brown', lwd = 2) # p DER
lines(d7_Si, si_one, col = 'blue') # Si NTE-also 1-ion DER
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 1) # I(d)

# The following chunk is just asking for trouble with unamed vectors
# Edward please fix
Prev=as.numeric(mix_data$Prev)-Y_0; SD=as.numeric(mix_data$SD)
errbar(100,Prev[1], yplus  = Prev[1] +  SD[1],
       yminus = Prev[1] -SD[1],
       pch = 19, cap = 0.05, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1)
#legend(x = "topleft", legend = "not 95%CI, SD", cex=0.3)

# ========================================================== #
#= Fig. 7B_final p-Si mix, TE, point & error, narrow ribbon =# 
# ========================================================== #
# We use the plot that takes adjustable parameter correlations into account
corr_fig_7B <- simulate_monte_carlo(n = 500, d7_new, LET_vals, ratios, model = "TE")
# The first argument, n, is the number of Monte Carlo repeats. Increase for
# greater accuracy. Decrease to speed up the program.
ci_data <- data.frame(dose = d7_new,
                      # Monte Carlo values
                      corrBottom = corr_fig_7B$monte_carlo[1, ],
                      corrTop = corr_fig_7B$monte_carlo[2, ], #
                      
                      # one-ion DERs 
                      si = calibrated_HZE_te_der(dose = d7_new, L = 70),
                      p = calibrated_low_LET_der(dose = d7_new, L = .4),
                      
                      # IEA baseline mixture DER I(d), denoted by id below
                      i = calculate_id(d7_new, LET_vals, ratios, model = "TE")[, 2])

plot(c(0, 100), c(0,.36), col = "white", bty = 'L', ann = FALSE) # Set plot area
polygon(x = c(d7_new, rev(d7_new)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", lwd = .4, border = "orange") # Narrow CI ribbon
si_te_one = calibrated_HZE_te_der(dose = d7_Si, L = 70)
p_one = calibrated_low_LET_der(dose = d7_p, L = .4)
lines(d7_p, p_one, col = 'brown', lwd = 2) # p DER
lines(d7_Si, si_te_one, col = 'blue') # Si TE-only 1-ion DER
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 1) # I(d)


Prev=as.numeric(mix_data$Prev)-Y_0; SD=as.numeric(mix_data$SD)
errbar(100,Prev[1], yplus  = Prev[1] +  SD[1],
       yminus = Prev[1] -SD[1],
       pch = 19, cap = 0.05, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1)

#==================================================================#
#== Fig8_final. 3-ion mix. Correlated vs Uncorr CI Overlay Plot ===#
#==================================================================#
# Consists of 3 HZE ions in our 5/20/2019 data set.
# Declare ratios and LET values for plot
ratios <- rep(1/3,3)
LET_vals <- c(70, 100, 193)
d8 <- c(0.1 * 0:9, 1:60)
d8t <- c(0.1 * 0:9, 1:20)
# We begin with the correlated plot
corr_fig_8 <- simulate_monte_carlo(n = 500, d8, LET_vals, ratios, model = "NTE")
# Comments for Fig. 10 apply with minor changes here and in some other lines
# We now calculate the uncorrelated Monte Carlo
uncorr_fig_8 <- simulate_monte_carlo(n = 500, d8, LET_vals, ratios, model = "NTE", vcov = FALSE)

ci_data <- data.frame(dose = d8,
                      # Monte Carlo values
                      corrBottom = corr_fig_8$monte_carlo[1, ],
                      corrTop = corr_fig_8$monte_carlo[2, ],
                      uncorrBottom = uncorr_fig_8$monte_carlo[1, ],
                      uncorrTop = uncorr_fig_8$monte_carlo[2, ],
                      
                      # DER values
                      si = calibrated_HZE_nte_der(dose = d8, L = 70),
                      ti = calibrated_HZE_nte_der(dose = d8, L = 100),
                      fe_six = calibrated_HZE_nte_der(dose = d8, L = 193),
                      i = calculate_id(d8, LET_vals, ratios, model = "NTE")[, 2])

# Plotting call. 
plot(c(0, 60), c(0, .450), col = "white", bty = 'l', ann = FALSE) # Next is broad CI ribbon
polygon(x = c(d8, rev(d8)), y = c(ci_data[, "uncorrTop"], rev(ci_data[, "uncorrBottom"])),
        xpd = -1, col = "aquamarine2", lwd = .5, border = "aquamarine2") # Wide CI

polygon(x = c(d8, rev(d8)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", border = "orange", lwd = .2) # Narrow CI

lines(d8t, calibrated_HZE_nte_der(dose = d8t, L = 70), col = 'black', lwd = 2)
lines(d8t, calibrated_HZE_nte_der(dose = d8t, L = 100), col = 'black', lwd = 2)
lines(d8t, calibrated_HZE_nte_der(dose = d8t, L = 193), col = 'black', lwd = 2)
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 3,lty=3) # I(d)
#errbar(49,.35, yplus  = .35 +  .05,
 #      yminus = .35 -.05,
  #     pch = 19, cap = 0.05, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1) 

#==================================================================#
#== Fig9_final. 8-ion mix. Correlated vs Uncorr CI Overlay Plot ===#
#==================================================================#
# Consists of all 8 HZE ions in our 5/20/2019 data set.
# Declare ratios and LET values for plot
ratios <- rep(1/8,8)
LET_vals <- c(20, 25, 70, 100, 193, 250, 464, 953)
d9 <- c(.001*0:9,0.1 * 1:9, 1:50)
d9t <- c(.001*0:9,0.1 * 1:9, 1:6,6.1,6.2,6.25)
d9B = c(.001*0:9,0.1 * 1:9, 1:10)
# We begin with the correlated plot
corr_fig_9 <- simulate_monte_carlo(n = 500, d9, LET_vals, ratios, model = "NTE")
# We now calculate the uncorrelated Monte Carlo
uncorr_fig_9 <- simulate_monte_carlo(n = 500, d9, LET_vals, ratios, model = "NTE", vcov = FALSE)

ci_data <- data.frame(dose = d9,
                      # Monte Carlo values
                      corrBottom = corr_fig_9$monte_carlo[1, ],
                      corrTop = corr_fig_9$monte_carlo[2, ],
                      uncorrBottom = uncorr_fig_9$monte_carlo[1, ],
                      uncorrTop = uncorr_fig_9$monte_carlo[2, ],
                      
                      # DER values
                      ne = calibrated_HZE_nte_der(dose = d9, L = 20),
                      o = calibrated_HZE_nte_der(dose = d9, L = 25),
                      si = calibrated_HZE_nte_der(dose = d9, L = 70),
                      ti = calibrated_HZE_nte_der(dose = d9, L = 100),
                      fe_six = calibrated_HZE_nte_der(dose = d9, L = 193),
                      fe_three = calibrated_HZE_nte_der(dose = d9, L = 250),
                      nb = calibrated_HZE_nte_der(dose = d9, L = 464),
                      la = calibrated_HZE_nte_der(dose = d9, L = 953),
                      i = calculate_id(d9, LET_vals, ratios, model = "NTE")[, 2])

# Plotting call. 
plot(c(0, 50), c(0, .40), col = "white", bty = 'l', ann = FALSE) # Next is broad CI ribbon
polygon(x = c(d9, rev(d9)), y = c(ci_data[, "uncorrTop"], rev(ci_data[, "uncorrBottom"])),
        xpd = -1, col = "aquamarine2", lwd = .5, border = "aquamarine2") # Wide CI

polygon(x = c(d9, rev(d9)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", border = "orange", lwd = .2) # Narrow CI

lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 20), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 25), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 70), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 100), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 193), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 250), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 464), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 953), col = 'blue', lwd = 1)
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 2) # I(d)
errbar(49,.35, yplus  = .35 +  .05,
       yminus = .35 -.05,
       pch = 19, cap = 0.05, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1) #OK up to here

#=====================================#
# panel B ============================#
#=====================================#
ratios <- rep(1/8,8)
LET_vals <- c(20, 25, 70, 100, 193, 250, 464, 953)
d9 <- c(.001*0:9,0.1 * 1:9, 1:50)
d9t <- c(.001*0:9,0.1 * 1:9, 1:6,6.1,6.2,6.25)
d9B = c(.001*0:9,0.1 * 1:9, 1:10)
# We begin with the correlated plot
corr_fig_9B <- simulate_monte_carlo(n = 500, d9B, LET_vals, ratios, model = "NTE")
# We now calculate the uncorrelated Monte Carlo
uncorr_fig_9B <- simulate_monte_carlo(n = 500, d9B, LET_vals, ratios, model = "NTE", vcov = FALSE)

ci_data <- data.frame(dose = d9B,
                      # Monte Carlo values
                      corrBottom = corr_fig_9B$monte_carlo[1, ],
                      corrTop = corr_fig_9B$monte_carlo[2, ],
                      uncorrBottom = uncorr_fig_9B$monte_carlo[1, ],
                      uncorrTop = uncorr_fig_9B$monte_carlo[2, ],
                      
                      # DER values
                      ne = calibrated_HZE_nte_der(dose = d9B, L = 20),
                      o = calibrated_HZE_nte_der(dose = d9B, L = 25),
                      si = calibrated_HZE_nte_der(dose = d9B, L = 70),
                      ti = calibrated_HZE_nte_der(dose = d9B, L = 100),
                      fe_six = calibrated_HZE_nte_der(dose = d9B, L = 193),
                      fe_three = calibrated_HZE_nte_der(dose = d9B, L = 250),
                      nb = calibrated_HZE_nte_der(dose = d9B, L = 464),
                      la = calibrated_HZE_nte_der(dose = d9B, L = 953),
                      i = calculate_id(d9B, LET_vals, ratios, model = "NTE")[, 2])

# Plotting call. 
plot(c(0, 10), c(0, .12), col = "white", bty = 'l', ann = FALSE) # Next is broad CI ribbon
polygon(x = c(d9B, rev(d9B)), y = c(ci_data[, "uncorrTop"], rev(ci_data[, "uncorrBottom"])),
        xpd = -1, col = "aquamarine2", lwd = .5, border = "aquamarine2") # Wide CI

polygon(x = c(d9B, rev(d9B)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", border = "orange", lwd = .2) # Narrow CI

lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 20), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 25), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 70), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 100), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 193), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 250), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 464), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 953), col = 'blue', lwd = 1)
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 2, lty = 2) # I(d)

#errbar(49,.35, yplus  = .35 +  .05,
       #yminus = .35 -.05,
       #pch = 19, cap = 0.05, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1)
#======================================================================#
#= Figqq_WebSup. Orders prevalences regardless of dose. Then compares =#
#= quartile cutpoint values with normal distribution for eyeball test =#
#======================obsolete =======================================#
# qq_in = c(low_LET_data[, "Prev"], HZE_data[,"Prev"]) # unordered
# qq_sequence = qqnorm(qq_in, ann = FALSE) # ordered and compared
# qqline(qq_in) # For Gaussian data all points would lie near this line, but
              # here: low tails are shorter, high dose tails longer
#==================================================================#
#== Fig9A_WebSup. 8-ion mix. Correlated ribbon; SEA ===#
#==================================================================#
# Consists of all 8 HZE ions in our 5/20/2019 data set.
# Declare ratios and LET values for plot
ratios <- rep(1/8,8)
LET_vals <- c(20, 25, 70, 100, 193, 250, 464, 953)
d9 <- c(.001*0:9,0.1 * 1:9, 1:50)
d9t <- c(.001*0:9,0.1 * 1:9, 1:6,6.1,6.2,6.25)
d9B = c(.001*0:9,0.1 * 1:9, 1:10)
# We use the correlated plot
corr_fig_9 <- simulate_monte_carlo(n = 500, d9, LET_vals, ratios, model = "NTE")

ci_data <- data.frame(dose = d9,
                      # Monte Carlo values
                      corrBottom = corr_fig_9$monte_carlo[1, ],
                      corrTop = corr_fig_9$monte_carlo[2, ],
                      
                      # DER values
                      o = calibrated_HZE_nte_der(dose = d9, L = 20),
                      ne = calibrated_HZE_nte_der(dose = d9, L = 25),
                      si = calibrated_HZE_nte_der(dose = d9, L = 70),
                      ti = calibrated_HZE_nte_der(dose = d9, L = 100),
                      fe_six = calibrated_HZE_nte_der(dose = d9, L = 193),
                      fe_three = calibrated_HZE_nte_der(dose = d9, L = 250),
                      nb = calibrated_HZE_nte_der(dose = d9, L = 464),
                      la = calibrated_HZE_nte_der(dose = d9, L = 953),
                      i = calculate_id(d9, LET_vals, ratios, model = "NTE")[, 2])

# Plotting call. 
plot(c(0, 50), c(0, .50), col = "white", bty = 'l', ann = FALSE) # Next is broad CI ribbon

polygon(x = c(d9, rev(d9)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", border = "orange", lwd = .2) # Narrow CI

lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 20), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 25), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 70), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 100), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 193), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 250), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 464), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 953), col = 'blue', lwd = 1)
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 2) # I(d)
errbar(49,.35, yplus  = .35 +  .05,
       yminus = .35 -.05,
       pch = 19, cap = 0.05, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1) #OK up to here
sea_der = calculate_SEA(d9, LET_vals, ratios)
lines(d9,sea_der,lty =2)