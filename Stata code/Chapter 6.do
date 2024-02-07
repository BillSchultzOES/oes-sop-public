
********************** 6.1

** Set root directory
cd "...\Stata code"

** Read in data for the fake experiment.
import delimited using "dat1_with_designs.csv", clear

egen ranky = rank(y), track
qui count
global ndat1 = r(N)

qui sum y1
gen y1s = y1 / r(sd)
qui sum y0
gen y0s = y0 / r(sd)

qui tabstat y1s y0s, stat(mean) save
global trueates = r(StatTotal)[1,2] - r(StatTotal)[1,1]
qui tabstat y1 y0, stat(mean) save
global trueate = r(StatTotal)[1,2] - r(StatTotal)[1,1]

**********

power twomeans 0, /// Requires specifying a control mean
n(1000) /// The total number of observations (default is equal division)
power(0.8) // The traditional power threshold of 80%

**********

qui power twomeans 0, n(1000) power(0.8) 
global delta = r(delta) // See: return list
di "$delta"

**********

power twomeans 0, diff($delta) power(0.8)

********************** 6.2

** OES Power Simulation Toolkit (Stata):
** 
** draw_from_design ---: Generate a simulated dataset (NOT RUN DIRECTLY)
** single_estimator ---: Draw data once and estimate results (NOT RUN DIRECTLY)
** replicate ----------: Repeat (generate -> estimate) many times
** evaluate_power -----: Evaluate power to detect non-zero effects.
** evaluate_mde -------: Find MDE, searching over range of effect sizes.
** evaluate_bias ------: Compute bias and other diagnostics.

***** DRAW FROM a hypothetical design *****

** Note: Must be modified by user
* Required set-up:
* - (1) Write code within this program to generate one draw of a simulated dataset.
*   - See our examples below, this is often simpler than it sounds!
*   - A built-in simulated treatment effect is not generally needed.
capture program drop draw_from_design
program define draw_from_design, nclass

	* Clear existing data
	clear
	
	** Replace the rest of the code inside this program with your own code
	
	* Sample size of 1000 observations
	set obs 1000
	
	* Generate simulated outcome
	gen y = rnormal(0, 1)
	
	* Generate simulated treatment (complete random assignment)
	qui count
	local ntreat = r(N)/2
	complete_ra x, m(`ntreat')

end

**** ESTIMATE results for a single simulated dataset ****

** Note: Must be modified by user
* Required set-up:
* - (1) Define the data generationprogram above.
* - (2) Write out the test you want to run in this program (using the simulated data).
capture program drop single_estimator
program define single_estimator, rclass

	* Check that design program exists
	quietly capture draw_from_design
	if _rc != 0 {
		di as error "Error: define data generation program (draw_from_design) first"
		exit
	}
	
	* Call the design program
	draw_from_design
	
	* Write out the desired estimation strategy
	reg y x, vce(hc2)

end

**** REPEAT (generation -> estimation) many times ****

** Note: Modification by user NOT NEEDED (just copy into your .do file)
* Required set-up:
* - (1) Define both programs above (data generation and a single_estimator)
* Inputs on use:
* - (1) Number of replicates (default = 200)
capture program drop replicate
program define replicate, rclass

	syntax[, reps(integer 200) ]

	* Check that design program exists
	quietly capture draw_from_design
	if _rc != 0 {
		di as error "Error: define data generation program (draw_from_design) first"
		exit
	}
	
	* Check that single_estimator program exists
	quietly capture single_estimator
	if _rc != 0 {
		di as error "Error: define estimation program (single_estimator) first"
		exit
	}

	* Save coefficients and SEs from each draw to memory.
	simulate ///
	_b _se, ///
	reps(`reps') nodots: ///
	single_estimator
	
	* Simulation indicator
	gen sim = _n
	
	* Modify var names of coefficients/SEs slightly
	foreach var of varlist _b* _se* {
		qui rename `var' `=substr("`var'", 2, .)'
	}	

	* Reshape to a format that makes the desired power calculation easier.
	qui reshape long b_ se_, i(sim) j(term) string
		
end

**** EVALUATE power of the design ****

** Note: Modification by user NOT NEEDED (just copy into your .do file)
* Required setup:
* - (1) Define all programs above
* - (2) Run replicate to get simulated coef/SE estimates in memory.
* Inputs on use:
* - (1) Hypothetical effects we want power estimates for (min, steps, and max)
*	- (Default: from 0 to 1 in steps of 0.01)
* - (2) Desired alpha (significance) level (default = 0.05)
capture program drop evaluate_power
program define evaluate_power, nclass

	syntax[, ///
	delta_min(real 0) ///
	delta_steps(real 0.01) ///
	delta_max(real 1) ///
	alpha(real 0.05) ]
	
	* Data to return to after each iteration
	tempfile restore_dat
	qui save `restore_dat', replace
	
	* Loop over specified effect sizes
	local i = 0
	forvalues n = `delta_min'(`delta_steps')`delta_max' {
		
		qui use `restore_dat', clear
		local ++i
		
		* Real and simulated statistics for each a given effect size
		gen delta = `n'
		gen real_t = b_/se_
		gen sim_t = (b_ + delta)/se_
		
		* Generate a p-value for each value of sim_t
		qui gen sim_p = .
		qui count
		forvalues v = 1/`r(N)' { // Loop over observations
			qui gen greaterequal = abs(real_t) >= abs(sim_t[`v'])
			qui sum greaterequal if term == term[`v'], meanonly
			qui replace sim_p = r(mean) if _n == `v'
			qui drop greaterequal
		}
		
		* Use these to get power for the given delta
		qui gen reject = sim_p <= `alpha'
		bysort term: egen power = mean(reject)
		collapse (mean) power, by(term delta)
		label var power ""
		
		* Save, and advance to the next delta
		if `i' == 1 {
			qui tempfile running_dat
			qui save `running_dat', replace
		}
		
		else {
			append using `running_dat'
			qui save `running_dat', replace
		}

	}
	
	* Open the result, replacing the data in memory
	qui use `running_dat', clear

end

**** EVALUATE the min. detectable effect ****

** Note: Modification by user NOT NEEDED (just copy into your .do file)
* Required setup:
* - (1) Define programs above
* - (2) Run replicate and evaluate_power
* Inputs on use:
* - (1) Minimum power desired
capture program drop evaluate_mde
program define evaluate_mde, nclass

	syntax[, min_power(real 0.8)]
	
	quietly {
		bysort term (power): gen above_min = power >= `min_power'
		drop if above_min == 0
		bysort term (delta): gen min = _n == 1
		drop if min == 0
		drop min above_min
	}

end

**** EVALUATE Bias for a particular term ****

** Note: Modification by user NOT NEEDED (just copy into your .do file)
* Required setup:
* - (1) Define programs above
* - (2) Run replicate
* Inputs on use:
* - (1) The name of the term to provide diagnosics for
* - (2) True parameter value (generally true ATE)
capture program drop evaluate_bias
program define evaluate_bias, nclass

	syntax, true_value(real) term(string) [restore_data]
	
	* Save data to return to in temporary file
	* (program includes option to turn this off)
	if "`restore_data'" != "" {
		tempfile restore_data
		qui save `restore_data', replace
	}

	* Subset to only the term in question
	qui keep if term == "`term'"
	
	* True parameter as variable
	qui gen true_value = `true_value'

	* Prepare variables to summarizetrue
	qui gen bias = b_ - true_value
	qui gen MSE = (b_ - true_value)^2
	qui gen conflow = b_ - (1.96 * se_) // normal approximation
	qui gen confhigh = b_ + (1.96 * se_) // normal approximation
	qui gen coverage = true_value >= conflow & true_value <= confhigh
	
	collapse ///
	(first) True = true_value ///
	(mean) Mean_Estimate = b_ Bias = bias MSE Coverage = coverage Mean_SE = se_ ///
	(sd) SD_Estimate = b_, ///
	by(term)

	list
	
	* Return to data?
	if "`restore_data'" != "" {
		qui use `restore_data', clear
	}
	
end

**********

** Analytical power estimates
* View power estimates for a range of effect sizes as a table
power twomeans 0, n(1000) diff(0.005(0.005)0.505)
* Or view as a graph instead
power twomeans 0, n(1000) diff(0.005(0.005)0.505) graph
* It's also possible to get the estimates as data in memory
* and write your own plotting code (e.g.: using twoway).
clear
svmat r(pss_table), names(col)
list in 1/5 // Illustration: power estimates are data in memory
keep diff power
rename diff delta
tempfile analytical
save `analytical', replace

** Computational power estimates, using the programs as defined above
replicate, reps(200) // Replicate data generation and estimation 200 times
evaluate_power, delta_min(0.005) delta_max(0.500) delta_steps(0.005) // Consider a range of effect sizes

** Merge computational with analytic estimates
rename power power_comp
keep if term == "x"
merge 1:1 delta using `analytical'
keep if _merge == 3
drop _merge

** Manual line plot
label var power "Analytical"
label var power_comp "Computational"
twoway ///
(line power delta) ///
(line power_comp delta), ///
legend(pos(6) rows(1)) ///
xtitle("Effect size") ytitle("Power")

**********

** 0. Define the programs used to simulate data and apply an estimator

* 0a: simulate data (update program below as needed)
* Output is a dataset in memory
capture program drop draw_from_design
program define draw_from_design, nclass

	* Clear existing data
	clear

	* Sample size of 1000 observations
	set obs 1000
	
	* Generate simulated outcome
	gen y = rnormal(0, 1)
	
	* Generate simulated treatment (complete random assignment)
	qui count
	local ntreat = r(N)/2
	complete_ra x, m(`ntreat')

end

* 0b: apply estimator (update program below as needed)
* Output is a dataset in memory and stored estimates
capture program drop single_estimator
program define single_estimator, rclass
	
	* Call the design program
	draw_from_design
	
	* Write out the desired estimation strategy
	reg y x, vce(hc2)

end

** 1/2. Replicate/Estimate:

* Output is a dataset of coefficients and SEs from each simulation.
replicate, reps(200) // Number of replications

** 3. Evaluate 

* Output is a set of power estimates in memory, one for each delta
evaluate_power, ///
delta_min(0.005) /// Smallest delta to consider
delta_max(0.500) /// Largest delta to consider
delta_steps(0.005) // Increments to apply

**********

** Basic elements of each simulated sample replicate
* Redefine data generation
capture program drop draw_from_design
program define draw_from_design, nclass

	* Clear existing data
	clear

	* Sample size of 1000 observations
	set obs 1000
	
	* Generate simulated outcome
	gen y = rbinomial(1, 0.25)
	
	* Generate simulated treatments
	qui count
	local ntreat = r(N)/2
	complete_ra x, m(`ntreat')
	complete_ra x2, m(`ntreat')

end

** Estimate main and interaction effects
* Redefine estimation
capture program drop single_estimator
program define single_estimator, rclass
	
	* Call the design program
	draw_from_design
	
	* Write out the desired estimation strategy
	* (note: the program currently does not correctly handle factor notation)
	gen x_int = x*x2
	reg y x x2 x_int, vce(hc2)

end

** Replicate estimates
replicate, reps(200)

** Evaluate power
evaluate_power, delta_min(0) delta_max(0.25) delta_steps(0.002)

**********

** Reshape to apply plotting code similar to above
reshape wide power, i(delta) j(term) string

** Line plot
label var powerx "x1"
label var powerx2 "x2"
label var powerx_int "x1 * x2"
twoway ///
(line powerx delta, lcolor(black) lpattern(solid)) ///
(line powerx2 delta, lcolor(black) lpattern(dash)) ///
(line powerx_int delta, lcolor(black) lpattern(dot)), ///
legend(pos(6) rows(1)) ///
xtitle("Effect size") ytitle("Power")

**********

** Basic elements of each simulated sample replicate
* Redefine data generation
capture program drop draw_from_design
program define draw_from_design, nclass

	* Clear existing data
	clear

	* Sample size of 1000 observations
	set obs 1000
	
	* Continuous covariate
	gen z1 = rnormal(0, 3)
	
	* Binary covariate
	gen z2 = rbinomial(1, 0.25)
	
	* Mean centered versions
	qui sum z1
	gen cz1 = z1 - r(mean)
	qui sum z2
	gen cz2 = z2 - r(mean)
	
	* Generate simulated treatment (10%)
	complete_ra x, prob(0.1)
	
	* Simulate y
	gen y = ((z1 + z2) * x) + (0.8 * z1) - (1 * z2) + rnormal()

end

**********

** With the lin estimator
capture program drop single_estimator
program define single_estimator, rclass
	
	* Call the design program
	draw_from_design
	
	* Write out the desired estimation strategy
	* (As noted above, this code isn't currently written for factor notation)
	gen x_cz2 = x * cz2
	gen x_cz1 = x * cz1
	reg y x x_cz2 x_cz1 cz1 cz2, vce(hc2)

end

* Replicate and save
replicate, reps(200)
evaluate_power, delta_min(0) delta_max(0.5) delta_steps(0.005)
keep if term == "x"
rename power power_lin
tempfile lin
save `lin', replace

** With no covariate adjustment
capture program drop single_estimator
program define single_estimator, rclass
	
	* Call the design program
	draw_from_design
	
	* Write out the desired estimation strategy
	reg y x, vce(hc2)

end

* Replicate and save
replicate, reps(200)
evaluate_power, delta_min(0) delta_max(0.5) delta_steps(0.005)
keep if term == "x"
rename power power_no
tempfile no
save `no', replace

** With linear, additive covariate adjustment
capture program drop single_estimator
program define single_estimator, rclass
	
	* Call the design program
	draw_from_design
	
	* Write out the desired estimation strategy
	reg y x cz1 cz2, vce(hc2)

end

* Replicate and merge
replicate, reps(200)
evaluate_power, delta_min(0) delta_max(0.5) delta_steps(0.005)
keep if term == "x"
merge 1:1 delta using `lin'
drop _merge
merge 1:1 delta using `no'
drop _merge

**********

** Line plot
label var power "Additive"
label var power_no "No covariates"
label var power_lin "Lin"
twoway ///
(line power delta, lcolor(black) lpattern(solid)) ///
(line power_lin delta, lcolor(black) lpattern(dash)) ///
(line power_no delta, lcolor(black) lpattern(dot)), ///
legend(pos(6) rows(1)) ///
xtitle("Effect size") ytitle("Power") ///
title("Power with Lin adjustment") yline(0.8)