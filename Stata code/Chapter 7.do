
**********************

** Set root directory
cd "...Stata code"

** Read in data for the fake experiment.
import delimited using "dat1_with_designs.csv", clear

**********************

** Create an ATE variable
gen tau = y1 - y0

** Plot the first 6 rows
list y1 y0 y z tau in 1/6

**********************

** Actual average treatment effect
di "Actual ATE"
qui sum tau, meanonly
di "`r(mean)'"

** Sample estimate of the ATE: manually
di "Sample estimate of ATE"
qui sum y if z==1, meanonly
local mean_z1 = r(mean)
qui sum y if z==0, meanonly
local mean_z0 = r(mean)
di `mean_z1'-`mean_z0'

** Sample estimate of the ATE: regression
*qui reg y z
*di _b[z]

**********************

** E.g., a Kolmogorov-Smirnov test for a
** difference in outcome distributions.
di "KS test for a difference in distributions"
ksmirnov y, by(z)

** E.g., a rank-based test for a "location shift"
** in the distribution of Y.
di "MWU rank-based test of a shift in distributions"
ranksum y, by(z)

**********************

** Baseline outcome rate.
** Assuming we might have control
** variables, let's compute this manually
** instead of using the regression intercept.
di "Average baseline outcome"
qui sum y if z == 0, meanonly
local mean_z0 = r(mean)
di "`mean_z0'"

** Sample ATE estimate
di "Sample ATE estimate"
qui reg y z
di _b[z]
local ate = _b[z]

** Percent change in baseline due to treatment
local p_chng = `ate'/`mean_z0'*100
di "Percent change = `p_chng'"

**********************

** Fit a regression model
qui reg y z
local df = e(df_r)

** Get estimates with HC2 errors instead of default SEs.
qui reg y z, vce(hc2)

** Compute t-stat from coef and SE
local t = _b[z]/_se[z]

** Compare to its null (of 0) sampling distribution.
** Compute two-sided p-value.
local p = (2 * ttail(`df', abs(`t')))

* View
di "`p'"

**********************

** Get CI based on HC2 SEs
qui reg y z, vce(hc2)
local df = e(df_r)
local ci_low = r(table)["ll","z"]
local ci_up = r(table)["ul","z"]

** Compute manually for comparison
local crit_val = invttail(`df', 0.025)
local ci_low2 = _b[z] - (_se[z]*`crit_val')
local ci_up2 = _b[z] + (_se[z]*`crit_val')

** View
di "`ci_low2' to `ci_up2'"

**********************

** Generate a treatment variable that is
** not associated with the outcome in any way.
set seed 1234
gen u = runiform()
gsort u
gen nullZ = 1 in 1/50
replace nullZ = 0 if missing(nullZ)

** Regress the outcome on this null treatment indicator.
qui reg y nullZ, vce(hc2)
local ci_low = r(table)["ll","nullZ"]
local ci_up = r(table)["ul","nullZ"]

** Get rough approximation of ex-post MDE.
** If, say, the SESI was 3, then this is
** greater than the MDE, which is good.
** But we still want to check the confidence interval.
local mde = _se[nullZ]*2.8
di "E.g.: SESI = 3"
di "MDE80 = `mde'"
di "CI = `ci_low' to `ci_up'"

**********************

** 90% CI
qui reg y nullZ, vce(hc2) level(90)
local df = e(df_r)
local ci_low = r(table)["ll","nullZ"]
local ci_up = r(table)["ul","nullZ"]

** Define equivalence region so that
** it is right on the edge of CI upper bound.
local eq_low = -2.25
local eq_high = 2.25
di "90% CI = `ci_low' to `ci_up'"
di "EQ region = `eq_low' to `eq_high'"

** First TOST one-sided p-value
local p1 = (_b[nullZ] - `eq_low')/_se[nullZ]
local p1 = ttail(`df', abs(`p1'))

** Second TOST one-sided p-value
local p2 = (_b[nullZ] - `eq_high')/_se[nullZ]
local p2 = 1 - ttail(`df', `p2')

** TOST overall p-value (95% confidence).
** Borderline, since 90% CI is almost at 
** edge of equivalence region!
local both_pvals `p1' `p2'
local max : subinstr local both_pvals " " ",", all
local max = max(`max')
di "TOST p-value = `max'"


