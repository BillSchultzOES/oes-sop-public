
**********************

** Set root directory
cd "...\Stata code"

**********************

** Read in data for the fake experiment.
import delimited using "dat1.csv", clear
rename v1 x

** Table of the first few observations.
list x y0 y1 z y id in 1/6, sep(6)

**********************

** y0 and y1 are the true underlying potential outcomes.
label var y0 "Control"
label var y1 "Treatment"
* ssc install stripplot
stripplot y0 y1, box vertical iqr whiskers(recast(rcap)) ytitle("Outcomes") variablelabels

**********************

qui tabstat y0 y1, stat(mean) save
* return list
* di r(StatTotal)
global trueATE = r(StatTotal)[1,2] - r(StatTotal)[1,1]

**********************

** Y is the observed outcome, Z is the observed treatment.
qui ttest y, by(z)
global estATE1 = round(r(mu_2) - r(mu_1), 0.001)
qui reg y z
global estATE2 = round(r(table)[1,1], 0.001)
di "Difference in means, $estATE1; OLS regression, $estATE2"
assert $estATE1 == $estATE2

**********************

** A program to re-assign treatment and recalculate the difference of means.
** Treatment was assigned without blocking or other structure, so we
** just permute or shuffle the existing treatment assignment vector.
capture program drop simEstAte
program define simEstAte, rclass sortpreserve

	version 18.0
	syntax varlist(min=1 max=1), control_outcome(varname) treat_outcome(varname)
	
	qui sum `varlist' // Get # treated units
	local numtreat = r(sum)
	
	tempvar rand // Randomly sort
	qui gen `rand' = runiform()
	sort `rand'
	
	tempvar Znew // Temporary var containing the new treatment
	qui gen `Znew' = 0
	qui replace `Znew' = 1 in 1/`numtreat'
	
	tempvar Ynew // Temporary var containing the new revealed outcome
	qui gen `Ynew' = (`Znew' * `treat_outcome') + ((1 - `Znew') * `control_outcome')
	
	qui ttest `Ynew', by(`Znew')
	return scalar estate = r(mu_2) - r(mu_1)
	
end

** Set up and perform the simulation
global sims 10000
set seed 12345
preserve
	qui simulate estate = r(estate), reps($sims): simEstAte z, control_outcome(y0) treat_outcome(y1)
	qui sum estate
	global seEstATEsim = `r(sd)'
restore

** The standard error of this estimate of the ATE (via simulation)
di "$seEstATEsim"

**********************

** True SE (Dunning Chap 6, Gerber and Green Chap 3, or Freedman, Pisani and Purves A-32).
** Requires knowing the true covariance between potential outcomes.
qui count
local N = r(N)
qui cor y0 y1, cov
local varc = r(C)[1,1]
local vart = r(C)[2,2]
local covtc = r(C)[1,2]
qui sum z
local nt = r(sum)
local nc = `N' - `nt'

** Gerber and Green, p.57, equation (3.4)
local varestATE = (((`varc' * `nt') / `nc') + ((`vart' * `nc') / `nt') + (2 * `covtc')) / (`N' - 1)
global seEstATETrue = sqrt(`varestATE')

**********************

** Feasible SE
qui sum y if z == 0
local varYc = r(sd) * r(sd)
qui sum y if z == 1
local varYt = r(sd) * r(sd)
local fvarestATE = (`N'/(`N'-1)) * ( (`varYt'/`nt') + (`varYc'/`nc') )
global estSEEstATE = sqrt(`fvarestATE')

**********************

** OLS SE
qui reg y z
global iidSE = _se[z] // Or: sqrt(e(V)["z","z"])

**********************

** Neyman SE (HC2)
qui reg y z, vce(hc2)
global NeymanSE = _se[z]  // Or: sqrt(e(V)["z","z"])

**********************

matrix compareSEs = J(1, 5, .)
matrix compareSEs[1, 1] = $seEstATEsim
matrix compareSEs[1, 2] = $estSEEstATE
matrix compareSEs[1, 3] = $seEstATETrue
matrix compareSEs[1, 4] = $iidSE
matrix compareSEs[1, 5] = $NeymanSE
matrix colnames compareSEs = "simSE" "feasibleSE" "trueSE" "olsIIDSE" "NeymanDesignSE"
matrix list compareSEs

**********************

** Define a function to calculate several SEs, given potential outcomes and treatment
capture program drop sePerfFn
program define sePerfFn, rclass sortpreserve

	version 18.0
	syntax varlist(min=1 max=1), control_outcome(varname) treat_outcome(varname)
	
	qui sum `varlist' // As in the program above
	local numtreat = r(sum)
	tempvar rand
	qui gen `rand' = runiform()
	sort `rand'
	tempvar Znew
	qui gen `Znew' = 0
	qui replace `Znew' = 1 in 1/`numtreat'
	tempvar Ynew
	qui gen `Ynew' = (`Znew' * `treat_outcome') + ((1 - `Znew') * `control_outcome')
	
	qui reg `Ynew' `Znew' // Regression now instead of ttest
	local iidSE = _se[`Znew']
	qui reg `Ynew' `Znew', vce(hc2)
	local NeymanSE = _se[`Znew']
	
	return scalar iidSE = `iidSE' // Prepare output
	return scalar NeymanSE = `NeymanSE'
	return scalar estATE = _b[`Znew']

end

** Perform a simulation using this function
set seed 12345

preserve
	qui simulate ///
	iidSE = r(iidSE) NeymanSE = r(NeymanSE) estATE = r(estATE), ///
	reps($sims): ///
	sePerfFn z, control_outcome(y0) treat_outcome(y1)
	qui sum estATE
	global simSE = r(sd)
	qui sum iidSE
	global estSEiid = r(mean)
	qui sum NeymanSE
	global estSENeyman = r(mean)
restore

di "Expected IID SE: $estSEiid"
di "Expected Neyman SE: $estSENeyman"
di "SIM SE: $simSE"
di "True SE: $seEstATETrue"

**********************

** Organize output
matrix compareCIs = J(2, 5, .)
matrix rownames compareCIs = "Approach 1 (diff. means)" "Approach 2 (OLS)"
matrix colnames compareCIs = "Est" "SE" "pvalue" "CI lower" "CI upper"

** A difference in means test assuming unequal variances
** (equivalent to the design-based estimator, as discussed above)
ttest y, by(z) unequal
local diffmeans = r(mu_2) - r(mu_1)
matrix compareCIs[1,1] = round(`diffmeans', 0.001)
matrix compareCIs[1,2] = round(r(se), 0.001)
matrix compareCIs[1,3] = round(r(p), 0.001)
matrix compareCIs[1,4] = round(`diffmeans' - (invttail(r(df_t), 0.025) * r(se)), 0.001)
matrix compareCIs[1,5] = round(`diffmeans' + (invttail(r(df_t), 0.025) * r(se)), 0.001)

** A regression-based approach (narrower intervals due to a different d.o.f)
reg y z, vce(hc2) // See: "matrix list r(table)"
matrix compareCIs[2,1] = round(r(table)[1, 1], 0.001)
matrix compareCIs[2,2] = round(r(table)[2, 1], 0.001)
matrix compareCIs[2,3] = round(r(table)[4, 1], 0.001)
matrix compareCIs[2,4] = round(r(table)[5, 1], 0.001)
matrix compareCIs[2,5] = round(r(table)[6, 1], 0.001)

matrix list compareCIs

**********************

** Define the grid to search through
numlist "-10(0.05)10" 
local grid `r(numlist)'
macro list _grid

** Simple RI program to use here
capture program drop ri_p
program define ri_p, rclass
	capture drop Zri
	qui complete_ra Zri, m(25)
	qui reg yg Zri
	return scalar taur = _b[Zri]
end

** Loop through values in this grid
tempfile realdat
save `realdat', replace
local i = 0
foreach g of local grid {
	
	local ++i
	
	** Create outcome for this g
	use `realdat', clear
	qui gen yg = y - `g' * z
	
	** Real regression for this g
	qui reg yg z
	local taug = _b[z]
	
	** RI inference for this g
	qui simulate taur = r(taur), ///
	reps(500) : ///
	ri_p
	replace taur = abs(taur) >= abs(`taug')
	collapse (mean) taur
	qui gen g = `g'
	
	** Save results
	if `i' == 1 {
		qui tempfile running
		qui save `running'
	}
	else {
		append using `running'
		qui save `running', replace
	}

	** Return to original data
	use `realdat', clear
	
}

** Load the results
use `running', clear

** Compute the CI
keep if taur > 0.05
qui sum g
global lower = r(min)
global upper = r(max)
di "$lower, $upper"