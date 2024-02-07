
********************** 5.1.1.1.1

** Set root directory
cd "...Stata code"

** Read in data for the fake experiment.
import delimited using "dat1_with_designs.csv", clear

**********

** Estimate the difference in means, assuming unequal variance
ttest y, by(z) unequal

**********

reg y z

********************** 5.1.1.1.2

** The permute command
set seed 12345

* Compare means
permute z z = _b[z], reps(1000) nodots: reg y z // OR: permtest2 y, by(z) simulate runs(1000)

* Rank test 1
egen ranky = rank(y)
permute z z = _b[z], reps(1000) nodots: reg ranky z

* Rank test 2 (note: exact, not approximate)
permute z z = r(z), reps(1000) nodots: ranksum y, by(z)
* OR: ranksum y, by(z) exact // exact, not approximate

**********

** The ritest command
* ssc install ritest
ritest z z = _b[z], nodots reps(1000): reg y z 

ritest z z = r(z), nodots reps(1000): ranksum y, by(z) 

**********

** Define a program to re-randomize a single time
** and return a desired test statistic.
capture program drop ri_draw
program define ri_draw, rclass

	** Using randomizr (Stata version)
	capture drop riZ
	qui sum z
	local zsum = r(sum)
	complete_ra riZ, m(`zsum')
	
	** Or manually	
		/*
		gen rand = runiform()
		sort rand
		qui sum z
		gen riZ = 1 in 1/r(sum)
		replace riZ = 0 if missing(riZ)
		drop rand
		*/
	
	** Return the test statistic of interest.
	** Simplest is the difference in means itself.
	qui reg y riZ
	return scalar riZ = _b[riZ]
	
end

** We'll use simulate to repeat this many times.
** Nested within preserve/restore to return to our data after.
preserve

	** Get the real test statistic
	qui reg y z
	local real_stat = _b[z]
	
	** Perform the simulation itself
	simulate ///
	riZ = r(riZ), ///
	reps(1000): ///
	ri_draw
	
	** Calculate the p-value.
	* How often do different possible treatment assignments,
	* under a sharp null, yield test statistics with a magnitude
	* at least as large as our real statistic?
	gen equal_or_greater = abs(riZ) >= abs(`real_stat')
	qui sum equal_or_greater
	local ri_p_manual = r(mean)
	
	** Compare to output from permute or ritest above
	di `ri_p_manual'

restore

********************** 5.1.1.2.1

** Make some binary outcomes
gen u = runiform()
gen v = runiform()
gen uv = u + v
gen y0bin = cond(u > 0.5, 1, 0) // control potential outcome
gen y1bin = cond(uv > 0.75, 1, 0) // treated potential outcome
gen ybin = (z * y1bin) + ((1 - z) * y0bin)
qui sum y1bin, meanonly
local y1mean = r(mean)
qui sum y0bin, meanonly
global truePropDiff = `y1mean' - r(mean)

** Estimate and view the difference in proportions
ttest ybin, by(z) unequal

**********

qui sum z
global nt = r(sum)
tempvar oneminus
gen `oneminus' = 1 - z
qui sum `oneminus'
global nc = r(sum)

** Find SE for difference of proportions.
qui ttest ybin, by(z)
local p1 = r(mu_2)
local p0 = r(mu_2)
local se1 = (`p1' * (1 - `p1'))/$nt
local se0 = (`p0' * (1 - `p0'))/$nc
global se_prop = round(sqrt(`se1' + `se0'), 0.0001)
local se_list $se_prop // Initialize a running list; used in the matrix below

** Find Neyman SE
qui sum ybin if z == 0
local varc_s = r(sd) * r(sd)
qui sum ybin if z == 1
local vart_s = r(sd) * r(sd)
global se_neyman = round(sqrt((`vart_s'/$nt) + (`varc_s'/$nc)), 0.0001)
local se_list `se_list' $se_neyman

** Find OLS SE
qui reg ybin z
global se_ols = round(_se[z], 0.0001) // See also: r(table) (return list) or e(V) (ereturn list)
local se_list `se_list' $se_ols

** Find Neyman SE (which are the HC2 SEs)
qui reg ybin z, vce(hc2)
global se_neyman2 = round(_se[z], 0.0001) 
qui ttest ybin, by(z) unequal // See: return list
global se_neyman3 = round(r(se), 0.0001)
local se_list `se_list' $se_neyman2 $se_neyman3

** Show SEs
matrix se_compare = J(1, 5, .)
matrix colnames se_compare = "diff in prop" "neyman1" "ols" "neyman2" "neyman3"
local i = 0
foreach l of local se_list {
	local ++i
	matrix se_compare[1, `i'] = `l'
}
matrix list se_compare

********************** 5.1.1.2.2

tabulate z ybin, exact

tabulate z ybin, chi2

* search emh
emh z ybin

**********

prtest ybin, by(z)

********************** 5.1.2

** Comparing only conditions 1 and 2 using ttest
ttest y if inlist(z4arms, "T1", "T2"), by(z4arms) unequal

**********

encode z4arms, gen(z4num)
reg y ib1.z4num, vce(hc2) // Set 1 as the reference category

********************** 5.1.2.1

** Get p-values but exclude intercept
* See also: search parmest
matrix pvals = r(table)["pvalue", 2..4] // save in a matrix
matrix pvalst = pvals' // transpose
svmat pvalst, names(col) // add matrix as data in memory

** Illustrate different corrections (or lack thereof)

* None
replace pvalue = round(pvalue, 0.0001)
list pvalue if !missing(pvalue)

* Bonferroni
* ssc install qqvalue
qqvalue pvalue if !missing(pvalue), method(bonferroni) qvalue(adj_p_bonf)
replace adj_p_bonf = round(adj_p_bonf, 0.0001)
list adj_p_bonf if !missing(pvalue)

* Holm
qqvalue pvalue if !missing(pvalue), method(holm) qvalue(adj_p_holm)
replace adj_p_holm = round(adj_p_holm, 0.0001)
list adj_p_holm if !missing(pvalue)

* Hochberg
qqvalue pvalue if !missing(pvalue), method(hochberg) qvalue(adj_p_hoch)
replace adj_p_hoch = round(adj_p_hoch, 0.0001)
list adj_p_hoch if !missing(pvalue) // FDR instead of FWER

**********

** First we need to estimate an anova model
anova y z4num

**********

* search tukeyhsd
* search qsturng
tukeyhsd z4num
* Or: pwcompare z4arms, mcompare(tukey) effects

********************** 5.2.2

** Preserve the running dataset used so far.
save main.dta, replace

** Keep a dataframe of select variables.
keep id y1 y0 cov*

** A dataset to represent a smaller experiment,
** or a cluster randomized experiment with few clusters
** (an experimental sample of 20 units).
** Only done for illustration here. In practice, we'll use the data
** generated in the R code for the sake of comparison.
set seed 12345
gen rand = runiform()
sort rand
keep in 1/20
drop rand

** Define a program to load one of these datasets
** and then randomly re-assign treatment.
capture program drop sample_from
program define sample_from, rclass

	syntax[ , smaller ///
		propsimtreat(real 0.5) ]
	
	* Which dataset to use?
	if "`smaller'" != "" import delimited using "smalldat1.csv", clear
	else import delimited using "popbigdat1.csv", clear
	
	* Make sure propsimtreat is a proportion
	if `propsimtreat' > 1 | `propsimtreat' < 0 {
		di as error "Check input: propsimtreat"
		exit
	}
	
	* Re-assign (simulated) treatment in this draw of the data
	complete_ra znew, prob(`propsimtreat')
	
	* No additional treatment effects
	* (assigning new potential outcomes)
	gen y_znew_0 = y0
	gen y_znew_1 = y1
	
	* Save the true ATE
	gen true_effect = y_znew_1 - y_znew_0
	qui sum true_effect, meanonly
	return scalar ATE = r(mean)
	
	* Get the revealed outcome
	gen ynew = (znew * y_znew_1) + ((1 - znew) * y_znew_0)

end

**********

** Define a program to apply various estimation strategies
** to a dataset drawn using the program defined above.
capture program drop apply_estimators
program define apply_estimators, rclass

	** Same arguments as above
	syntax[, smaller ///
		propsimtreat(real 0.5) ]
		
	** Call the program above
	sample_from, `smaller' propsimtreat(`propsimtreat')
	return scalar ATE = r(ATE)

	** CovAdj0: Lm, No Covariates
	qui reg ynew znew, vce(hc2)
	return scalar CovAdj0_est = _b[znew]
	return scalar CovAdj0_p = r(table)["pvalue", "znew"]
	
	** CovAdj1: Lm, Correct Covariate
	qui reg ynew znew cov2, vce(hc2)
	return scalar CovAdj1_est = _b[znew]
	return scalar CovAdj1_p = r(table)["pvalue", "znew"]
	
	** CovAdj2: Lm, Mixed Covariates
	qui reg ynew znew cov1-cov8, vce(hc2)
	return scalar CovAdj2_est = _b[znew]
	return scalar CovAdj2_p = r(table)["pvalue", "znew"]
	
	** CovAdj3: Lm, Wrong Covariates
	qui reg ynew znew cov1 cov3-cov6, vce(hc2)
	return scalar CovAdj3_est = _b[znew]
	return scalar CovAdj3_p = r(table)["pvalue", "znew"]
	
	** CovAdj4: Lin, Mixed Covariates
	qui reg ynew znew cov1-cov8 // to make a sample indicator
	gen samp = e(sample) // ensure correct obs are used in mean-centering
	foreach var of varlist cov1-cov8 {
		qui sum `var' if samp == 1, meanonly
		qui gen mc_`var' = `var' - `r(mean)' if samp == 1
	}
	qui reg ynew i.znew##c.(mc_*), vce(hc2)
	return scalar CovAdj4_est = _b[1.znew]
	return scalar CovAdj4_p = r(table)["pvalue", "1.znew"]
	drop mc_* samp
	
	** CovAdj5: Lin, Correct Covariate
	qui reg ynew znew cov2
	gen samp = e(sample)
	foreach var of varlist cov2 {
		qui sum `var' if samp == 1, meanonly
		qui gen mc_`var' = `var' - `r(mean)' if samp == 1
	}
	qui reg ynew i.znew##c.(mc_*), vce(hc2)
	return scalar CovAdj5_est = _b[1.znew]
	return scalar CovAdj5_p = r(table)["pvalue", "1.znew"]
	drop mc_* samp

end

**********

** Call the program once just to make a list of scalars to save
apply_estimators, smaller
local rscalars: r(scalars) // Save all scalars in r() to a local
local to_store "" // Update them in a loop to work properly in simulate
foreach item of local rscalars {
  local to_store "`to_store' `item' = r(`item')"
}

** Summarize characteristics of the smaller-sample designs
set seed 12345
simulate ///
`to_store', ///
reps(200): ///
apply_estimators, smaller

** Create summary matrix
qui des, short // saves number of columns to r()
matrix diagnosands = J((r(k) - 1)/2, 6, .)
matrix rownames diagnosands = "Lm, No Covariates" "Lm, Correct Covariate" "Lm, Mixed Covariates" "Lm, Wrong Covariates" "Lin, Mixed Covariates" "Lin, Correct Covariate"
matrix colnames diagnosands = "Mean estimand" "Mean estimate" "Bias" "SD Estimate" "RMSE" "Power"

** Calculate quantities to include
** (https://declaredesign.org/r/declaredesign/reference/declare_diagnosands.html)
local row = 0
forvalues i = 0/5 {
	
	local ++row
	
	* Estimand
	qui sum ATE, meanonly
	matrix diagnosands[`row', 1] = r(mean)
	
	* Estimate
	qui sum CovAdj`i'_est, meanonly
	matrix diagnosands[`row', 2] = r(mean)
	
	* Bias
	qui gen biascalc = CovAdj`i'_est - ATE
	qui sum biascalc, meanonly
	matrix diagnosands[`row', 3] = r(mean)
	drop biascalc
	
	* SD estimate (based on population variance, no n-1 in the variance denom.)
	qui sum CovAdj`i'_est
	qui gen sdcalc = (CovAdj`i'_est - r(mean))^2
	qui sum sdcalc
	matrix diagnosands[`row', 4] = sqrt(r(sum)/r(N))
	drop sdcalc
	
	* RMSE
	gen atediff = (CovAdj`i'_est - ATE)^2
	qui sum atediff, meanonly
	matrix diagnosands[`row', 5] = sqrt(r(mean))
	drop atediff
	
	* Power
	gen rejectnull = CovAdj`i'_p <= 0.05
	qui sum rejectnull, meanonly
	matrix diagnosands[`row', 6] = r(mean)
	drop rejectnull
	
}

* View the results
matrix list diagnosands

**********

** Summarize characteristics of the large-sample designs
apply_estimators
local rscalars : r(scalars)
local to_store ""
foreach item of local rscalars {
  local to_store "`to_store' `item' = r(`item')"
}
simulate ///
`to_store', ///
reps(200): ///
apply_estimators

** Create summary matrix
qui des, short // get number of columns in r()
matrix diagnosands = J((r(k) - 1)/2, 6, .)
matrix rownames diagnosands = "Lm, No Covariates" "Lm, Correct Covariate" "Lm, Mixed Covariates" "Lm, Wrong Covariates" "Lin, Mixed Covariates" "Lin, Correct Covariate"
matrix colnames diagnosands = "Mean estimand" "Mean estimate" "Bias" "SD Estimate" "RMSE" "Power"

** Calculate quantities to include
local row = 0
forvalues i = 0/5 {
	
	local ++row
	
	* Estimand
	qui sum ATE, meanonly
	matrix diagnosands[`row', 1] = r(mean)
	
	* Estimate
	qui sum CovAdj`i'_est, meanonly
	matrix diagnosands[`row', 2] = r(mean)
	
	* Bias
	qui gen biascalc = CovAdj`i'_est - ATE
	qui sum biascalc, meanonly
	matrix diagnosands[`row', 3] = r(mean)
	drop biascalc
	
	* SD estimate (based on population variance, no n-1 in the variance denom.)
	qui sum CovAdj`i'_est
	qui gen sdcalc = (CovAdj`i'_est - r(mean))^2
	qui sum sdcalc
	matrix diagnosands[`row', 4] = sqrt(r(sum)/r(N))
	drop sdcalc
	
	* RMSE
	gen atediff = (CovAdj`i'_est - ATE)^2
	qui sum atediff, meanonly
	matrix diagnosands[`row', 5] = sqrt(r(mean))
	drop atediff
	
	* Power
	gen rejectnull = CovAdj`i'_p <= 0.05
	qui sum rejectnull, meanonly
	matrix diagnosands[`row', 6] = r(mean)
	drop rejectnull
	
}

* View the results
matrix list diagnosands

********************** 5.2.3

** Define a program to perform Rosenbaum estimation
* (with all variables specified as in regress).
capture program drop rosenbaum_est
program define rosenbaum_est, rclass

	* Assumed: outcome treatment list_of_covariates.
	syntax varlist
	* Alternative 1: syntax, outcome(varname) treatment(varname) covar(varlist)
	* Alternative 2: syntax varlist, outcome(varname) treatment(varname)
	
	* Pull out first var of varlist as outcome
	tokenize `varlist'
	local outcome `1'
	macro shift
	
	* Pull out second var of varlist as treatment
	local treat_cov `*'
	tokenize `treat_cov'
	local treatment `1'
	macro shift
	
	* The rest are treated as covariates
	local covar `*'
	
	* Estimation
	reg `outcome' `covar'
	tempvar yhat
	predict `yhat', xb
	tempvar resid
	gen `resid' = `outcome' - `yhat'
	reg `resid' `treatment', vce(hc2)
	
	* Output
	return scalar out_est = _b[`treatment']
	return scalar out_p = r(table)["pvalue", "`treatment'"]

end

** Program to apply this estimator to the generated data:
capture program drop apply_estimators_2
program define apply_estimators_2, rclass

	syntax[, smaller ///
		propsimtreat(real 0.5) ]
		
	** Call the sample generation / random assignment program
	sample_from, `smaller' propsimtreat(`propsimtreat')
	return scalar ATE = r(ATE)

	** CovAdj6: Resid, Correct
	qui rosenbaum_est ynew znew cov2
	return scalar CovAdj6_est = r(out_est)
	return scalar CovAdj6_p = r(out_p)
	
	** CovAdj7: Resid, Mixed
	qui rosenbaum_est ynew znew cov1-cov8
	return scalar CovAdj7_est = r(out_est)
	return scalar CovAdj7_p = r(out_p)
	
	** CovAdj8: Resid, Incorrect
	qui rosenbaum_est ynew znew cov1 cov3-cov6
	return scalar CovAdj8_est = r(out_est)
	return scalar CovAdj8_p = r(out_p)

end

**********

** Summarize characteristics of the smaller-sample designs
apply_estimators_2, smaller // Looking only at the rosenbaum estimators now
local rscalars : r(scalars)
local to_store ""
foreach item of local rscalars {
  local to_store "`to_store' `item' = r(`item')"
}
simulate ///
`to_store', ///
reps(200): ///
apply_estimators_2, smaller

** Create summary matrix
qui des, short
matrix diagnosands = J((`r(k)' - 1)/2, 6, .)
matrix rownames diagnosands = "Resid, Correct" "Resid, Mixed" "Resid, Incorrect"
matrix colnames diagnosands = "Mean estimand" "Mean estimate" "Bias" "SD Estimate" "RMSE" "Power"

** Calculate quantities to include
local row = 0
forvalues i = 6/8 {
	
	local ++row
	
	* Estimand
	qui sum ATE, meanonly
	matrix diagnosands[`row', 1] = r(mean)
	
	* Estimate
	qui sum CovAdj`i'_est, meanonly
	matrix diagnosands[`row', 2] = r(mean)
	
	* Bias
	qui gen biascalc = CovAdj`i'_est - ATE
	qui sum biascalc, meanonly
	matrix diagnosands[`row', 3] = r(mean)
	drop biascalc
	
	* SD estimate (based on population variance, no n-1 in the variance denom.)
	qui sum CovAdj`i'_est
	qui gen sdcalc = (CovAdj`i'_est - r(mean))^2
	qui sum sdcalc
	matrix diagnosands[`row', 4] = sqrt(r(sum)/r(N))
	drop sdcalc
	
	* RMSE
	gen atediff = (CovAdj`i'_est - ATE)^2
	qui sum atediff, meanonly
	matrix diagnosands[`row', 5] = sqrt(r(mean))
	drop atediff
	
	* Power
	gen rejectnull = CovAdj`i'_p <= 0.05
	qui sum rejectnull, meanonly
	matrix diagnosands[`row', 6] = r(mean)
	drop rejectnull
	
}

* View the results
matrix list diagnosands

**********

** Summarize characteristics of the larger-sample designs
apply_estimators_2 // Looking only at the rosenbaum estimators now
local rscalars : r(scalars)
local to_store ""
foreach item of local rscalars {
  local to_store "`to_store' `item' = r(`item')"
}
di "`to_store'"
simulate ///
`to_store', ///
reps(200): ///
apply_estimators_2

** Create summary matrix
qui des, short
matrix diagnosands = J((`r(k)' - 1)/2, 6, .)
matrix rownames diagnosands = "Resid, Correct" "Resid, Mixed" "Resid, Incorrect"
matrix colnames diagnosands = "Mean estimand" "Mean estimate" "Bias" "SD Estimate" "RMSE" "Power"

** Calculate quantities to include
local row = 0
forvalues i = 6/8 {
	
	local ++row
	
	* Estimand
	qui sum ATE, meanonly
	matrix diagnosands[`row', 1] = r(mean)
	
	* Estimate
	qui sum CovAdj`i'_est, meanonly
	matrix diagnosands[`row', 2] = r(mean)
	
	* Bias
	qui gen biascalc = CovAdj`i'_est - ATE
	qui sum biascalc, meanonly
	matrix diagnosands[`row', 3] = r(mean)
	drop biascalc
	
	* SD estimate (based on population variance, no n-1 in the variance denom.)
	qui sum CovAdj`i'_est
	qui gen sdcalc = (CovAdj`i'_est - r(mean))^2
	qui sum sdcalc
	matrix diagnosands[`row', 4] = sqrt(r(sum)/r(N))
	drop sdcalc
	
	* RMSE
	gen atediff = (CovAdj`i'_est - ATE)^2
	qui sum atediff, meanonly
	matrix diagnosands[`row', 5] = sqrt(r(mean))
	drop atediff
	
	* Power
	gen rejectnull = CovAdj`i'_p <= 0.05
	qui sum rejectnull, meanonly
	matrix diagnosands[`row', 6] = r(mean)
	drop rejectnull
	
}

* View the results
matrix list diagnosands

* Return to main data from before the covariance adjustment simulations
use main.dta, clear
erase main.dta

********************** 5.4.2

** Data simulation provided for illustration.
** In practice, will use data generated by the R code,
** for the sake of comparison.

** Create block sizes and create block weights
clear
local sizes 8 20 30 40 50 60 70 80 100 800 // Block sample sizes
set obs 10 // Number of blocks
qui gen bf = . // Categorical block var
qui gen nb = . // Number in block
local i = 0
foreach size of local sizes {
	local ++i
	replace nb = `size' if _n == `i'
	replace bf = `i' if _n == `i'
}
gen lambda = runiform(1, 2000) // For poisson generation below
expand nb
sort bf

** x1 is a covariate that strongly predicts the outcome without treatment
set seed 2201
gen x1 = rpoisson(lambda)

** The treatment effect varies by block size (sqrt(nb) because nb has such a large range.)
bysort bf: egen sd_x1 = sd(x1) 
gen y0 = (sd_x1 * x1) + rchi2(1)
bysort bf: egen p5 = pctile(y0), p(5)
replace y0 = 0 if y0 <= p5
bysort bf: egen sd_y0 = sd(y0)
gen tauib = -(sd_y0)*sqrt(nb) + rnormal(0, sd_y0)
gen y1 = y0 + tauib
replace y1 = 0 if y1 <= 0
qui reg y0 i.bf
di round(`e(r2)', 0.0001)
drop sd_y0 sd_x1 p5

**********

** Define a program to load base dataset and randomly re-assign treatment
capture program drop sample_from
program define sample_from, rclass
	
	* Load simulated data from R
	import delimited using "blocksimdat.csv", clear
	
	* Get # treated for each block
	levelsof nb, local(nblevels)
	local numtreat
	local i = 0
	foreach l of local nblevels {
		
		local ++i // block indices (1, 2, 3,... 10)
		
		* even block indices (2nd, 4th, etc.) get 10% treated
		if mod(`i',2) == 0 {
			local temp = `l'/10 
			local numtreat `numtreat' `temp'
		}
		
		* odd get 50% treated
		else {
			local temp = `l'/2
			local numtreat `numtreat' `temp'
		}
		
	}
	
	* Re-assign (simulated) treatment in this draw of the data
	capture drop __0*
	block_ra znew, block_var(bf) block_m(`numtreat') replace
	
	* No additional treatment effects
	* (assigning new potential outcomes)
	gen y_znew_0 = y0
	gen y_znew_1 = y1

	* Get true ATE in r()
	gen true_effect = y_znew_1 - y_znew_0
	sum true_effect, meanonly
	return scalar ATE = r(mean)
	
	* Revealed outcome
	gen ynew = (znew * y_znew_1) + ((1 - znew) * y_znew_0)
	
	* Now add individual-level weights to the data
	* (under complete_ra within blocks, these will
	* be the same across re-randomizations).
	bysort bf: egen pib = mean(znew) // Prob. of treatment assignment
	bysort bf: egen nTb = total(znew) // Number treated
	gen nCb = nb - nTb // Number control
	gen nbwt = (znew/pib) + ((1-znew)/(1-pib)) // Unbiased regression weight
	gen hbwt = nbwt * (pib * (1 - pib)) // Precision regression weight

end

** Draw data using this program once and save the true ATE
set seed 2201
sample_from
global trueATE1 = round(r(ATE), 0.01)

** For illustration: how to prepare block-level data as in the R code.
** In practice, we'll use data from the R code for comparison.

** Prepare to create a block level dataset, with block level weights.
bysort bf znew: egen treatmean = mean(ynew) if znew == 1 // Treat and control means
bysort bf znew: egen controlmean = mean(ynew) if znew == 0
bysort bf znew: egen treatsd = sd(ynew) if znew == 1 // Treat and control SD
bysort bf znew: egen controlsd = sd(ynew) if znew == 0
qui count
gen total_n = r(N)

** Collapse to block-level
collapse ///
(mean) treatmean controlmean treatsd controlsd y1 y0 pb = znew ///
(first) nb nTb nCb total_n, ///
by(bf)

** Additional variable preparation
gen taub = treatmean - controlmean // Observed effect (variance below)
gen truetaub = y1 - y0 // True effect
gen treatvar = treatsd * treatsd
gen controlvar = controlsd * controlsd
gen estvartaub = (nb/(nb-1)) * (treatvar/nTb) + (controlvar/nCb)
gen nbwt = nb / total_n // Block size weight
gen pbwt = pb * (1 - pb) // pb is proportion treated
gen hbwt2 = nbwt * pbwt // Precision weights
gen hbwt3 = pbwt * nb
gen hbwt = (2*(nCb * nTb)/(nTb + nCb))
qui sum hbwt3
gen greenlabrule = 20 * hbwt3 / r(sum)

** All of these expressions of the harmonic mean weight are the same
qui sum hbwt
gen double hbwt01 = round(hbwt / r(sum), 0.000001) // Precision issues
qui sum hbwt2
gen double hbwt02 = round(hbwt2 / r(sum), 0.000001)
qui sum hbwt3
gen double hbwt03 = round(hbwt3 / r(sum), 0.000001)
assert hbwt01 == hbwt02
assert hbwt01 == hbwt03

** What is the "true" ATE?
gen truetaubweighted = truetaub * nbwt
qui sum truetaubweighted
global trueATE2 = round(r(sum), 0.01)
assert $trueATE1 == $trueATE2 // After rounding
** We could define the following as an estimand, too.
* gen truetaubweighted2 = truetaub * hbwt01
* qui sum truetaubweighted2
* global trueATE3 = round(`r(sum)', 0.01)

**********

** Load block-level data from the R code
import delimited using "datB.csv", clear

** simple_block
gen taubweighted = taub * nbwt
qui sum taubweighted
global ate_nbwt1 = r(sum)
tempvar calc_se
gen `calc_se' = nbwt^2 * estvartaub
qui sum `calc_se'
global ate_nbwt1se = sqrt(r(sum))

** design-based diffmeans estimator (return to observation level data)
import delimited using "blocksimdat2.csv", clear // dat2 in the R code
rename (y z ntb ncb) (ynew znew nTb nCb) // Adjust names to match illustrative code above
egen blockmean = mean(ynew), by(bf znew)
egen blocksd = sd(ynew), by(bf znew)
bysort bf (znew): gen diffmean = blockmean[_N] - blockmean[1] // with sorting, treat - control
bysort bf (znew): gen variance = ((blocksd[_N]^2)/nTb[_N]) + ((blocksd[1]^2)/nCb[1])
qui sum diffmean, meanonly
global ate_nbwt2 = r(mean)
qui count
gen forblockedse = (nb/r(N))^2 * variance
bysort bf: replace forblockedse = . if _n != 1
qui sum forblockedse
global ate_nbwt2se = sqrt(r(sum))

** lmlinbyhand
tabulate bf, gen(block_)
drop block_1
qui reg ynew znew i.bf
gen samp = e(sample) // Ensure correct sample used
foreach var of varlist block_* {
	qui sum `var' if samp == 1, meanonly
	gen mc_`var' = `var' - `r(mean)'
}
qui reg ynew i.znew##c.(mc_block_*), vce(hc2)
global ate_nbwt3 = _b[1.znew]
global ate_nbwt3se = _se[1.znew]

** regwts
qui reg ynew znew [aw = nbwt], vce(hc2)
global ate_nbwt6 = _b[znew]
global ate_nbwt6se = _se[znew]

** List all
di "simple_block = $ate_nbwt1"
di "diffmeans = $ate_nbwt2"
di "lmlinbyhand = $ate_nbwt3"
di "regwts = $ate_nbwt6"

**********

** Comparing the Standard Errors
di "simple_block = $ate_nbwt1se"
di "diffmeans = $ate_nbwt2se"
di "lmlinbyhand = $ate_nbwt3se"
di "regwts = $ate_nbwt6se"

**********

** simple_block (return to block level data)
import delimited using "datB.csv", clear
capture drop taubweighted_prec
gen taubweighted_prec = taub * hbwt01
qui sum taubweighted_prec
global ate_hbwt1 = r(sum)
capture drop __0*
tempvar calc_se2
gen `calc_se2' = hbwt01^2 * estvartaub
qui sum `calc_se2'
global ate_hbwt1se = sqrt(r(sum))

** lm_fixed_effects1 (return to observation level data)
import delimited using "blocksimdat2.csv", clear // dat2 in the R code
rename (y z ntb ncb) (ynew znew nTb nCb) // Adjust names to match illustrative code above
qui reg ynew znew i.bf, vce(hc2)
global ate_hbwt2 = _b[znew]
global ate_hbwt2se = _se[znew]

** lm_fixed_effects2 (SEs differ from R variant of this)
qui areg ynew znew, absorb(bf) vce(hc2)
global ate_hbwt3 = _b[znew]
global ate_hbwt3se = _se[znew]

** direct_wts
qui reg ynew znew [aw = hbwt], vce(hc2)
global ate_hbwt4 = _b[znew]
global ate_hbwt4se = _se[znew]

** demeaned
bysort bf: egen meanz = mean(znew)
bysort bf: egen meany = mean(ynew)
gen demeany = ynew - meany
gen demeanz = znew - meanz
qui reg demeany demeanz, vce(hc2)
global ate_hbwt5 = _b[demeanz]
global ate_hbwt5se = _se[demeanz]

** List all
di "simple_block = $ate_hbwt1"
di "lm_fixed_effects1 = $ate_hbwt2"
di "lm_fixed_effects2 = $ate_hbwt3"
di "direct_wts = $ate_hbwt4"
di "demeaned = $ate_hbwt5"

**********

** Comparing the Standard Errors
di "simple_block = $ate_hbwt1se"
di "lm_fixed_effects1 = $ate_hbwt2se"
di "lm_fixed_effects2 = $ate_hbwt3se"
di "direct_wts = $ate_hbwt4se"
di "demeaned = $ate_hbwt5se"

**********

** Define a program to apply various estimation strategies
** to a dataset drawn via sample_from.
capture program drop apply_estimators
program define apply_estimators, rclass
		
	** Call the data generation program defined above
	sample_from
	return scalar ATE = r(ATE)
	
	** E0: Ignores Blocks, OLS SE
	qui reg ynew znew
	return scalar estnowtIID_est = _b[znew]
	return scalar estnowtIID_se = _se[znew]
	return scalar estnowtIID_p = r(table)["pvalue", "znew"]
	
	** E1: Ignores Blocks, Design (HC2) SE
	qui reg ynew znew, vce(hc2)
	return scalar estnowtHC2_est = _b[znew]
	return scalar estnowtHC2_se = _se[znew]
	return scalar estnowtHC2_p = r(table)["pvalue", "znew"]
	
	** E2: Design-based difference in means, Design SE
	egen blockmean = mean(ynew), by(bf znew)
	egen blocksd = sd(ynew), by(bf znew)
	bysort bf (znew): gen diffmean = blockmean[_N] - blockmean[1] // with sorting, treat - control
	bysort bf (znew): gen variance = ((blocksd[_N]^2)/nTb[_N]) + ((blocksd[1]^2)/nTb[1])
	qui sum diffmean, meanonly
	local est = r(mean)
	return scalar estnbwt1_est = `est'
	qui count
	local N = r(N)
	gen forblockedse = (nb/`N')^2 * variance
	bysort bf: replace forblockedse = . if _n != 1
	qui sum forblockedse
	local se = sqrt(r(sum))
	return scalar estnbwt1_se = `se'
	qui levelsof bf, local(b)
	local J: word count `b'
	local df = `N' - (2*`J')
	return scalar estnbwt1_p = (2 * ttail(`df', abs(`est'/`se')))
	
	** E3: Treatment Interaction with Block Indicators, Design SE
	qui tabulate bf, gen(block_)
	drop block_1
	qui reg ynew znew i.bf
	gen samp = e(sample)
	foreach var of varlist block_* {
		qui sum `var' if samp == 1, meanonly
		gen mc_`var' = `var' - r(mean)
	}
	qui reg ynew i.znew##c.(mc_block_*), vce(hc2)
	return scalar estnbwt2_est = _b[1.znew]
	return scalar estnbwt2_se = _se[1.znew]
	return scalar estnbwt2_p = r(table)["pvalue", "1.znew"]
	
	** E5: Least Squares with Block Size Weights, Design SE
	qui reg ynew znew [aw = nbwt], vce(hc2)
	return scalar estnbwt4_est = _b[znew]
	return scalar estnbwt4_se = _se[znew]
	return scalar estnbwt4_p = r(table)["pvalue", "znew"]

	** E6: Precision Weights via Fixed Effects, Design SE
	qui reg ynew znew i.bf, vce(hc2)
	return scalar esthbwt1_est = _b[znew]
	return scalar esthbwt1_se = _se[znew]
	return scalar esthbwt1_p = r(table)["pvalue", "znew"]
	
	** E7: Precision Weights via Demeaning, Design SE
	qui areg ynew znew, absorb(bf) vce(hc2)
	return scalar esthbwt2_est = _b[znew]
	return scalar esthbwt2_se = _se[znew]
	return scalar esthbwt2_p = r(table)["pvalue", "znew"]
	
	** E8: Direct Precision Weights, Design SE
	qui reg ynew znew [aw = hbwt], vce(hc2)
	return scalar esthbwt3_est = _b[znew]
	return scalar esthbwt3_se = _se[znew]
	return scalar esthbwt3_p = r(table)["pvalue", "znew"]
	
	** E9: Direct Demeaning, Design SE
	bysort bf: egen meanz = mean(znew)
	bysort bf: egen meany = mean(ynew)
	gen demeany = ynew - meany
	gen demeanz = znew - meanz
	qui reg demeany demeanz, vce(hc2)
	return scalar esthbwt4_est = _b[demeanz]
	return scalar esthbwt4_se = _se[demeanz]
	return scalar esthbwt4_p = r(table)["pvalue", "demeanz"]
	
end

** Perform the simulation
apply_estimators
local rscalars : r(scalars)
local to_store ""
foreach item of local rscalars {
  local to_store "`to_store' `item' = r(`item')"
}
simulate ///
`to_store', ///
reps(200): ///
apply_estimators

** Create summary matrix
qui des, short
di (`r(k)' - 1)/3
matrix diagnosands = J((`r(k)' - 1)/3, 6, .)
matrix rownames diagnosands = "E0: Ignores Blocks, OLS SE" "E1: Ignores Blocks, Design (HC2) SE" "E2: Diff Means BW, Design SE" "E3: Lin BW, Design SE" "E5: Direct BW, Design SE" "E6: FE PW, Design SE" "E7: Demean PW, Design SE" "E8: Direct PW, Design SE" "E9: Direct Demeaning, Design SE"
matrix colnames diagnosands = "Mean estimand" "Mean estimate" "Bias" "SD Estimate" "MeanSE" "Coverage"

** Get variable name roots to loop through
local roots estnowtIID estnowtHC2 estnbwt1 estnbwt2 estnbwt4 esthbwt1 esthbwt2 esthbwt3 esthbwt4
 
** Calculate quantities to include
local row = 0
foreach l of local roots {
	
	local ++row
	
	* Estimand
	qui sum ATE, meanonly
	matrix diagnosands[`row', 1] = `r(mean)'
	
	* Estimate
	qui sum `l'_est, meanonly
	matrix diagnosands[`row', 2] = `r(mean)'
	
	* Bias
	qui gen biascalc = `l'_est - ATE
	qui sum biascalc, meanonly
	matrix diagnosands[`row', 3] = r(mean)
	drop biascalc
	
	* SD estimate (based on population variance, no n-1 in the variance denom.)
	qui sum `l'_est
	qui gen sdcalc = (`l'_est - r(mean))^2
	qui sum sdcalc
	matrix diagnosands[`row', 4] = sqrt(r(sum)/r(N))
	drop sdcalc
	
	* Mean SE
	qui sum `l'_se
	matrix diagnosands[`row', 5] = `r(mean)'
	
	* Coverage (z-stat CIs to limit estimates stored above)
	qui gen conflow = `l'_est - (1.96 * `l'_se)
	qui gen confhigh = `l'_est + (1.96 * `l'_se)
	qui gen inrange = ATE <= confhigh & ATE >= conflow
	qui sum inrange
	matrix diagnosands[`row', 6] = r(mean)
	drop conf* inrange
	
}

* View the results
matrix list diagnosands

**********

** Perform this simulation again with rank-transformed potential outcomes

** Redefine sampling program to consider ranked potential outcomes
capture program drop sample_from
program define sample_from, rclass
	
	* Load simulated data from R
	import delimited using "blocksimdat.csv", clear
	
	* Get # treated for each block
	levelsof nb, local(nblevels)
	local numtreat
	local i = 0
	foreach l of local nblevels {
		
		local ++i // block indices (1, 2, 3,... 10)
		
		* even block indices (2nd, 4th, etc.) get 10% treated
		if mod(`i',2) == 0 {
			local temp = `l'/10 
			local numtreat `numtreat' `temp'
		}
		
		* odd get 50% treated
		else {
			local temp = `l'/2
			local numtreat `numtreat' `temp'
		}
		
	}
	
	* Re-assign (simulated) treatment in this draw of the data
	capture drop __0*
	block_ra znew, block_var(bf) block_m(`numtreat') replace
	
	* No additional treatment effects, but rank POs
	egen y_znew_0 = rank(y0), unique // Different than how ties are handled in R.
	egen y_znew_1 = rank(y1), unique // R rank command averages by default.
	
	* Get true ATE in r()
	gen true_effect = y_znew_1 - y_znew_0
	sum true_effect, meanonly
	return scalar ATE = r(mean)
	
	* Revealed outcome
	gen ynew = (znew * y_znew_1) + ((1 - znew) * y_znew_0)
	
	* Now add individual-level weights to the data
	* (under complete_ra within blocks, these will
	* be the same across re-randomizations).
	bysort bf: egen pib = mean(znew) // Prob. of treatment assignment
	bysort bf: egen nTb = total(znew) // Number treated
	gen nCb = nb - nTb // Number control
	gen nbwt = (znew/pib) + ((1-znew)/(1-pib)) // Unbiased regression weight
	gen hbwt = nbwt * (pib * (1 - pib)) // Precision regression weight

end

** Estimators program above can then be re-used
apply_estimators
local rscalars : r(scalars)
local to_store ""
foreach item of local rscalars {
  local to_store "`to_store' `item' = r(`item')"
}
simulate ///
`to_store', ///
reps(200): ///
apply_estimators

** Create summary matrix
qui des, short
di (`r(k)' - 1)/3
matrix diagnosands = J((`r(k)' - 1)/3, 6, .)
matrix rownames diagnosands = "E0: Ignores Blocks, OLS SE" "E1: Ignores Blocks, Design (HC2) SE" "E2: Diff Means BW, Design SE" "E3: Lin BW, Design SE" "E5: Direct BW, Design SE" "E6: FE PW, Design SE" "E7: Demean PW, Design SE" "E8: Direct PW, Design SE" "E9: Direct Demeaning, Design SE"
matrix colnames diagnosands = "Mean estimand" "Mean estimate" "Bias" "SD Estimate" "MeanSE" "Coverage"

** Get variable name roots to loop through
local roots estnowtIID estnowtHC2 estnbwt1 estnbwt2 estnbwt4 esthbwt1 esthbwt2 esthbwt3 esthbwt4
 
** Calculate quantities to include
local row = 0
foreach l of local roots {
	
	local ++row
	
	* Estimand
	qui sum ATE, meanonly
	matrix diagnosands[`row', 1] = `r(mean)'
	
	* Estimate
	qui sum `l'_est, meanonly
	matrix diagnosands[`row', 2] = `r(mean)'
	
	* Bias
	qui gen biascalc = `l'_est - ATE
	qui sum biascalc, meanonly
	matrix diagnosands[`row', 3] = r(mean)
	drop biascalc
	
	* SD estimate (based on population variance, no n-1 in the variance denom.)
	qui sum `l'_est
	qui gen sdcalc = (`l'_est - r(mean))^2
	qui sum sdcalc
	matrix diagnosands[`row', 4] = sqrt(r(sum)/r(N))
	drop sdcalc
	
	* Mean SE
	qui sum `l'_se
	matrix diagnosands[`row', 5] = `r(mean)'
	
	* Coverage (z-stat CIs to limit estimates stored above)
	qui gen conflow = `l'_est - (1.96 * `l'_se)
	qui gen confhigh = `l'_est + (1.96 * `l'_se)
	qui gen inrange = ATE <= confhigh & ATE >= conflow
	qui sum inrange
	matrix diagnosands[`row', 6] = r(mean)
	drop conf* inrange
	
}

* View the results
matrix list diagnosands

**********

** Perform this simulation again with equal assignment probabilities.

** Redefine sampling program once more.
capture program drop sample_from
program define sample_from, rclass
	
	* Load simulated data from R
	import delimited using "blocksimdat.csv", clear

	* Re-assign (simulated) treatment in this draw of the data.
	capture drop __0*
	block_ra znew, block_var(bf) replace
	
	* No additional treatment effects, but rank POs
	egen y_znew_0 = rank(y0), unique // Different than how ties are handled in R.
	egen y_znew_1 = rank(y1), unique // R rank command averages by default.
	
	* Get true ATE in r()
	gen true_effect = y_znew_1 - y_znew_0
	sum true_effect, meanonly
	return scalar ATE = r(mean)
	
	* Revealed outcome
	gen ynew = (znew * y_znew_1) + ((1 - znew) * y_znew_0)
	
	* Now add individual-level weights to the data
	* (under complete_ra within blocks, these will
	* be the same across re-randomizations).
	bysort bf: egen pib = mean(znew) // Prob. of treatment assignment
	bysort bf: egen nTb = total(znew) // Number treated
	gen nCb = nb - nTb // Number control
	gen nbwt = (znew/pib) + ((1-znew)/(1-pib)) // Unbiased regression weight
	gen hbwt = nbwt * (pib * (1 - pib)) // Precision regression weight

end

** Estimators program above can then be re-used
apply_estimators
local rscalars : r(scalars)
local to_store ""
foreach item of local rscalars {
  local to_store "`to_store' `item' = r(`item')"
}
simulate ///
`to_store', ///
reps(200): ///
apply_estimators

** Create summary matrix
qui des, short
di (`r(k)' - 1)/3
matrix diagnosands = J((`r(k)' - 1)/3, 6, .)
matrix rownames diagnosands = "E0: Ignores Blocks, OLS SE" "E1: Ignores Blocks, Design (HC2) SE" "E2: Diff Means BW, Design SE" "E3: Lin BW, Design SE" "E5: Direct BW, Design SE" "E6: FE PW, Design SE" "E7: Demean PW, Design SE" "E8: Direct PW, Design SE" "E9: Direct Demeaning, Design SE"
matrix colnames diagnosands = "Mean estimand" "Mean estimate" "Bias" "SD Estimate" "MeanSE" "Coverage"

** Get variable name roots to loop through
local roots estnowtIID estnowtHC2 estnbwt1 estnbwt2 estnbwt4 esthbwt1 esthbwt2 esthbwt3 esthbwt4
 
** Calculate quantities to include
local row = 0
foreach l of local roots {
	
	local ++row
	
	* Estimand
	qui sum ATE, meanonly
	matrix diagnosands[`row', 1] = `r(mean)'
	
	* Estimate
	qui sum `l'_est, meanonly
	matrix diagnosands[`row', 2] = `r(mean)'
	
	* Bias
	qui gen biascalc = `l'_est - ATE
	qui sum biascalc, meanonly
	matrix diagnosands[`row', 3] = r(mean)
	drop biascalc
	
	* SD estimate (based on population variance, no n-1 in the variance denom.)
	qui sum `l'_est
	qui gen sdcalc = (`l'_est - r(mean))^2
	qui sum sdcalc
	matrix diagnosands[`row', 4] = sqrt(r(sum)/r(N))
	drop sdcalc
	
	* Mean SE
	qui sum `l'_se
	matrix diagnosands[`row', 5] = `r(mean)'
	
	* Coverage (z-stat CIs to limit estimates stored above)
	qui gen conflow = `l'_est - (1.96 * `l'_se)
	qui gen confhigh = `l'_est + (1.96 * `l'_se)
	qui gen inrange = ATE <= confhigh & ATE >= conflow
	qui sum inrange
	matrix diagnosands[`row', 6] = r(mean)
	drop conf* inrange
	
}

* View the results
matrix list diagnosands

********************** 5.5

** Modify dat2 to use for our clustering examples
import delimited using "blocksimdat2.csv", clear // dat2 in the R code
gen cluster = b

** Randomly assign half of the clusters to treatment and half to control
set seed 12345
cluster_ra zcluster, cluster_var(cluster)
table zcluster cluster

**********

** Update the sampling program for a basic clustered design
capture program drop sample_from
program define sample_from, rclass
	
	* Load baseline data
	import delimited using "blocksimdat2.csv", clear // dat2 in the R code
	gen cluster = b
	
	* Re-assign (simulated) treatment in this draw of the data
	cluster_ra zcluster, cluster_var(cluster)
	
	* No additional treatment effects
	gen y_zcluster_0 = y0
	gen y_zcluster_1 = y1
	
	* Get true ATE in r()
	gen true_effect = y_zcluster_1 - y_zcluster_0
	sum true_effect, meanonly
	return scalar ATE = r(mean)
	
	* Revealed outcome
	gen ynew = (zcluster * y_zcluster_1) + ((1 - zcluster) * y_zcluster_0)
	
	** Get weights for one estimation strategy used below
	qui loneway ynew cluster
	gen varweight = 1 / ( r(sd_b)^2 + ( r(sd_w)^2 / nb ) ) 

end

**********

** Get the treatment assignment used in the R code
import delimited using "datCluster.csv", clear

** Stata's default clustered errors
reg y zcluster, vce(cluster cluster)
local CR1 = _se[zcluster]

** CR2 errors via reg_sandwich
* ssc install reg_sandwich
reg_sandwich y zcluster, cluster(cluster)
local CR2 = _se[zcluster]

** Stata's default is CR1, not CR2
macro list _CR1 _CR2
assert `CR1' != `CR2'

********************** 5.5.1

** Define estimators that can be repeated in the simulation below
capture program drop apply_estimators
program define apply_estimators, rclass

	syntax[ , datCluster]
	
	if "`datCluster'" == "" {
		** Call the data generation program defined above
		sample_from
		return scalar ATE = r(ATE)
	}
	
	else {
		** Get the treatment assignment used in the R code
		import delimited using "datCluster.csv", clear
		rename y ynew
	}
	
	** C0: Ignores Clusters, IID SE
	qui reg ynew zcluster
	return scalar estC0_est = _b[zcluster]
	return scalar estC0_se = _se[zcluster]
	return scalar estC0_p = r(table)["pvalue", "zcluster"]
	
	** C1: Ignores Clusters, HC2 SE
	qui reg ynew zcluster, vce(hc2)
	return scalar estC1_est = _b[zcluster]
	return scalar estC1_se = _se[zcluster]
	return scalar estC1_p = r(table)["pvalue", "zcluster"]
	
	** C2: OLS, CR1 SE (Stata)
	qui reg ynew zcluster, vce(cluster cluster)
	return scalar estC2_est = _b[zcluster]
	return scalar estC2_se = _se[zcluster]
	return scalar estC2_p = r(table)["pvalue", "zcluster"]

	** C3: OLS, CR2 SE
	qui reg_sandwich ynew zcluster, cluster(cluster)
	return scalar estC3_est = _b[zcluster]
	return scalar estC3_se = _se[zcluster]
	local p = (2 * ttail(e(dfs)[1,1], abs(_b[zcluster]/_se[zcluster])))
	return scalar estC3_p = `p'

	** C6: OLS with Cluster Size Weights, CR2 SE
	qui reg_sandwich ynew zcluster [pw = nb], cluster(cluster)
	return scalar estC6_est = _b[zcluster]
	return scalar estC6_se = _se[zcluster]
	local p = (2 * ttail(e(dfs)[1,1], abs(_b[zcluster]/_se[zcluster])))
	return scalar estC6_p = `p'
	
	** C7: OLS with Cluster Control and Var. Weights, CR2 SE
	qui reg_sandwich ynew zcluster nb [pw = varweight], cluster(cluster)
	return scalar estC7_est = _b[zcluster]
	return scalar estC7_se = _se[zcluster]
	local p = (2 * ttail(e(dfs)[1,1], abs(_b[zcluster]/_se[zcluster])))
	return scalar estC7_p = `p'
	
	** C8: OLS Clusters with adj for cluster size, CR2 RCSE
	qui reg_sandwich ynew zcluster nb, cluster(cluster)
	return scalar estC8_est = _b[zcluster]
	return scalar estC8_se = _se[zcluster]
	local p = (2 * ttail(e(dfs)[1,1], abs(_b[zcluster]/_se[zcluster])))
	return scalar estC8_p = `p'
	
	** C9: OLS Clusters with adj for cluster size, CR2 RCSE
	qui reg ynew zcluster nb
	gen sample = e(sample)
	qui sum nb if sample == 1
	gen mc_nb = nb - r(mean) if sample == 1
	gen interaction = zcluster * mc_nb
	qui reg_sandwich ynew zcluster mc_nb interaction, cluster(cluster)
	return scalar estC9_est = _b[zcluster]
	return scalar estC9_se = _se[zcluster]
	local p = (2 * ttail(e(dfs)[1,1], abs(_b[zcluster]/_se[zcluster])))
	return scalar estC9_p = `p'
	
end

** Apply each of these estimators to the data and review the results.
* Built in option to fit all the models to datCluster (from the R code)
* for illustration rather than a data draw with a random treatment variable 
apply_estimators, datCluster

**********

** Simulate the performance of these estimators
apply_estimators
local rscalars : r(scalars)
local to_store ""
foreach item of local rscalars {
  local to_store "`to_store' `item' = r(`item')"
}
simulate ///
`to_store', ///
reps(200): ///
apply_estimators

** Create summary matrix
qui des, short
di (`r(k)' - 1)/3
matrix diagnosands = J((`r(k)' - 1)/3, 6, .)
matrix rownames diagnosands = "C0: IID SE" "C1: HC2 SE" "C2: CR1 SE" "C3: CR2 SE" "C6: Cluster Size Weights" "C7: Var Weights" "C8: Cluster size control" "C9: Lin cluster size adjust"
matrix colnames diagnosands = "Mean estimand" "Mean estimate" "Bias" "SD Estimate" "Power" "Coverage"

** Get variable name roots to loop through
local roots estC0 estC1 estC2 estC3 estC6 estC7 estC8 estC9
 
** Calculate quantities to include
local row = 0
foreach l of local roots {
	
	local ++row
	
	* Estimand
	qui sum ATE, meanonly
	matrix diagnosands[`row', 1] = `r(mean)'
	
	* Estimate
	qui sum `l'_est, meanonly
	matrix diagnosands[`row', 2] = `r(mean)'
	
	* Bias
	qui gen biascalc = `l'_est - ATE
	qui sum biascalc, meanonly
	matrix diagnosands[`row', 3] = r(mean)
	drop biascalc
	
	* SD estimate (based on population variance, no n-1 in the variance denom.)
	qui sum `l'_est
	qui gen sdcalc = (`l'_est - r(mean))^2
	qui sum sdcalc
	matrix diagnosands[`row', 4] = sqrt(r(sum)/r(N))
	drop sdcalc
	
	* Power
	gen rejectnull = `l'_p <= 0.05
	qui sum rejectnull, meanonly
	matrix diagnosands[`row', 5] = r(mean)
	drop rejectnull
	
	* Coverage (z-stat CIs to limit estimates stored above)
	qui gen conflow = `l'_est - (1.96 * `l'_se)
	qui gen confhigh = `l'_est + (1.96 * `l'_se)
	qui gen inrange = ATE <= confhigh & ATE >= conflow
	qui sum inrange
	matrix diagnosands[`row', 6] = r(mean)
	drop conf* inrange
	
}

* View the results
matrix list diagnosands

********************** 5.5.2

** Return to data from Chapter 4
import delimited using "dat1_with_designs.csv", clear

** Break any relationship between treatment and outcomes by permuting
** or shuffling the treatment variable. This means that H0, the null,
** of no effects is true.
capture program drop checkFP
program define checkFP, rclass

	** Specify if you want CR2, rather than CR1 (Stata)
	syntax[, se_CR2]
	capture drop newZ
	cluster_ra newZ, cluster_var(buildingid)
	if "`se_CR2'" == "" {
		reg_sandwich y newZ, cluster(buildingid)
		local p = (2 * ttail(e(dfs)[1,1], abs(_b[newZ]/_se[newZ])))
		return scalar nullp = `p'
	}
	else {
		reg y newZ, vce(cluster buildingid)
		local p = (2 * ttail(e(df_r), abs(_b[newZ]/_se[newZ])))
		return scalar nullp = `p'
	}

end

** Apply the function above 1000 times (CR2)
preserve
	set seed 123
	simulate ///
	nullp = r(nullp), ///
	reps(1000) nodots: ///
	checkFP, se_CR2
	
	qui gen reject = nullp <= 0.05
	qui sum reject
	di "fprateCR205 `r(mean)'"
restore

** Apply the function above 1000 times (Stata)
preserve
	set seed 123
	simulate ///
	nullp = r(nullp), ///
	reps(1000) nodots: ///
	checkFP
	
	qui gen reject = nullp <= 0.05
	qui sum reject
	di "fprateStata05 `r(mean)'"
restore
