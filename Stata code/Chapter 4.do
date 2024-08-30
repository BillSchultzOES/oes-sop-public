
********************** 4.1

** Set root directory
*cd "...Stata code"
cd "C:\Users\WilliamBSchultz\Documents\oes-sop-public\Stata code"

**********

** Start with a small experiment with only 10 units
clear
global n = 10
set obs $n

** Set a random seed for replicability
set seed 12345

** Simulation using functions from the randomizr package
* ssc install randomizr
simple_ra trt_coinflip
* Or, e.g.: gen trt_coinflip = rbinomial(1, 0.5)
complete_ra trt_urn
/* Or, e.g.:
local num_treat = $n/2
gen rand = runiform()
sort rand
gen trt_urn = 1 in 1/`num_treat'
replace trt_urn = 0 if missing(trt_urn)
drop rand
*/

** Add more informative labels
label define tc 0 "C" 1 "T"
label values trt_coinflip tc
label values trt_urn tc

** Coin flipping does not guarantee half and half treated and control.
** Drawing from an urn, guarantees half treated and control.
table trt_coinflip
table trt_urn

********************** 4.2

** Read in data for the fake experiment.
import delimited using "dat1.csv", clear

qui count
global N = r(N)
set seed 12345

** Two equal arms
complete_ra z2armEqual
label define tc 0 "C" 1 "T"
label values z2armEqual tc

** Two unequal arms: .25 chance of treatment (.75 chance of control)
complete_ra z2armUnequalA, prob(0.25)
label values z2armUnequalA tc
qui sum z2armUnequalA
global expected = $N/4
assert r(sum) == $expected
complete_ra z2armUnequalB, m($expected)
label values z2armUnequalB tc

** Four equal arms
local count_list : di _dup(4) "$expected " // List of sample sizes for each group
macro list _count_list
complete_ra z4arms, m_each(`count_list')

table z2armEqual
table z2armUnequalA
table z2armUnequalB
table z4arms

********************** 4.3

** Two equal arms, adding a second cross treatment
complete_ra z2armEqual2
label var z2armEqual2 "Treatment 2"
label var z2armEqual "Treatment 1"
label val z2armEqual2 tc
table z2armEqual z2armEqual2

********************** 4.?

local Anum = $N/2
gen blockID = .
tempvar rand
gen `rand' = runiform()
sort `rand'
replace blockID = 1 in 1/`Anum'
replace blockID = 2 if missing(blockID)
label define blocklab 1 "Block A" 2 "Block B"
label values blockID blocklab
table blockID

**********

block_ra z2armBlocked, block_var(blockID) replace
label val z2armBlocked tc
table blockID z2armBlocked
capture drop __00* // Clean up

**********

block_ra z2armBlockedUneqProb, block_var(blockID) block_prob(0.5 0.25) replace
label val z2armBlockedUneqProb tc
table blockID z2armBlockedUneqProb
capture drop __00* // Clean up

** Unit Testing
qui count if z2armBlockedUneqProb == 1 & blockID == 2
global NumTreatedB = `r(N)'
qui count if blockID == 2
global ExpectedNumTreatedB = `r(N)'/4
assert ($NumTreatedB == ceil($ExpectedNumTreatedB)) | ($NumTreatedB == floor($ExpectedNumTreatedB))

********************** 4.4.1

** For example, make three groups from the cov2 variable
egen cov2cat = cut(cov2), group(3) label
table cov2cat
bysort cov2cat : sum cov2 // Note: divides into intervals differently from R

** And we can make blocks that are the same on two covariates
qui sum cov1, d
gen cov1bin = cond(cov1 > `r(p50)', 1, 0) // Similar to ifelse() in R
decode cov2cat, generate(string_cov2cat)
gen blockV2 = string(cov1bin) + " " + string_cov2cat
table blockV2

** And then assign within these blocks
set seed 12345
block_ra zblockV2, block_var(blockV2)
label val zblockV2 tc
table blockV2 zblockV2
capture drop __00* // Clean up

********************** 4.5

** Make an indicator for cluster membership
qui count
global ndat1 = r(N)
egen buildingID = seq(), from(1) to(10) block(1)
set seed 12345
cluster_ra zcluster, cluster_var(buildingID)
table zcluster buildingID

********************** 4.6a

** load data to match R code
import delimited "ch4finaldat1.csv", clear

** List covariates to evaluate
global covlist cov1 cov2

** Define convenience program to perform calculations.
** May be applied to the whole dataset, cluster-level aggregates,
** or individual blocks. Assumes only 2 arms. See use below.
capture program drop bstats_2arm
program define bstats_2arm, rclass sortpreserve byable(onecall)

	* varname = evaluate balance for what covariate?
	* tr = treatment var name
	* n = total sample size; only use when applied to a subset (Blocks)
	syntax varname [if] [in], tr(varname) ///
		[ n(string) psd(string) ttestopt(string) ]
	
	* Get list of by vars, passed from "by varlist:" or "bysort varlist:"
	local by "`_byvars'"
	
	* Identify the correct sample to use (if/in)
	marksample touse
	
	* Further by var setup.
	* Confirm something was passed.
	capture confirm variable `by'
	
	* If not, treat all obs as in same group.
	if _rc != 0 {
		tempvar group
		qui gen `group' = 1 if `touse'	
	}
	
	* If something was, set this up as a grouping var.
	else {		
		tempvar group
		qui egen `group' = group(`by') if `touse'		
	}
	
	* In either case, get the levels as a macro.
	qui levelsof `group' if `touse', local(by_levels)
	local len : word count `by_levels'
	
	* Make rownames, used below
	local rownames
	foreach l of local by_levels {
		local rownames `rownames' "`l'"
	}
	
	* Get n as sample size in memory, if not specified
	if "`n'"=="" {
		qui count if `touse'
		local n = r(N)
	}
	
	* Get psd as pooled sd of varname, if not specified.
	* Though varname indicated above, still mapped to `varlist'.
	if "`psd'"=="" {
		qui sum `varlist' if `touse', d
		local psd = r(sd)
	}

	* Loop through those levels, calculating
	* desired stats and saving in a matrix.
	matrix stats = J(`len', 5, .)
	matrix rownames stats = `rownames'
	matrix colnames stats = "dm" "smd" "vr" "prop" "group_n"
	local i = 0
	foreach l of local by_levels {
		
		* Update index
		local ++i
		
		* Difference in means (treatment level 2 is greater value)
		qui ttest `varlist' if `touse' & `group' == `l', by(`tr') `ttestopt'
		local dm = r(mu_2) - r(mu_1)
		matrix stats[`i', 1] = `dm'
		
		* Standardize using provided SD
		local smd = `dm'/`psd'
		matrix stats[`i', 2] = `smd'
		
		* Variance ratio
		qui ttest `varlist' if `touse' & `group' == `l', by(`tr') `ttestopt'
		local vr = (r(sd_2)^2)/(r(sd_1)^2)
		matrix stats[`i', 3] = `vr'
		
		* Proportion of total sample in subset
		* (1 if n was not set)
		qui count if `touse' & `group' == `l'
		local prop = r(N)/`n'
		local group_n = r(N)
		matrix stats[`i', 4] = `prop'
		matrix stats[`i', 5] = `group_n'
			
	}
	
	* Prepare output
	return matrix stats = stats

end

** Iterate through this covariate list, building a table
local clen : word count $covlist
matrix btab = J(`clen', 6, .)
local rownames
foreach g of global covlist {
	local rownames `rownames' "`g'"
}
matrix rownames btab = `rownames'
matrix colnames btab = "smd" "vratio" "cluster_smd" "cluster_vratio" "block_smd" "block_vratio"
local k = 0
foreach var of varlist $covlist {
	
	* approach 1: overall mean differences and SMDs
	local ++k
	qui sum `var', d
	local pooled = r(sd)
	bstats_2arm `var', tr(z2armequal) psd(`pooled')
	matrix btab[`k', 1] = r(stats)[1,"smd"]
	matrix btab[`k', 2] = r(stats)[1,"vr"]
	
	* approach 2: cluster-level mean differences and SMD.
	* approach 1 might be applied for clustered designs as well.
	preserve
		qui collapse (mean) cov1 cov2 (first) zcluster, by(buildingid)
		bstats_2arm `var', tr(zcluster) psd(`pooled') ttestopt("reverse")
	restore
	matrix btab[`k', 3] = r(stats)[1,"smd"]
	matrix btab[`k', 4] = r(stats)[1,"vr"]
	
	* approach 3: weighted avg of in-block mean differences and SMDs.
	preserve
		bysort blockv2: bstats_2arm `var', tr(zblockv2) psd(`pooled') n(100)
		qui clear
		qui svmat r(stats), names(col)
		qui sum smd [iw=prop]
		local blocked_smd = r(mean)
		qui sum vr [iw=prop]
		local blocked_vr = r(mean)
	restore
	matrix btab[`k', 5] = `blocked_smd'
	matrix btab[`k', 6] = `blocked_vr'
	
}

** View output
matrix list btab

********************** 4.6b

** Prepare data for this exercise, save prior data
tempfile restore
save `restore', replace
local treatlist z2armequal zcluster zblockv2
foreach l of local treatlist {
	replace `l' = cond(`l' == "T", "1", "0")
	destring `l', replace
}

** Complete RA omnibus balance check
qui reg z2armequal cov1 cov2, vce(hc2)
test cov1=cov2=0
local waldp = r(p)
local waldf = r(F)

** Repeat using randomization inference
capture program drop ri_draw
program define ri_draw, rclass
	capture drop riZ
	complete_ra riZ
	qui reg riZ cov1 cov2, vce(hc2)
	qui test cov1=cov2=0
	return scalar ri_F = r(F)
end
preserve
	simulate ///
	ri_F = r(ri_F), ///
	reps(1000) nodots: ///
	ri_draw
	gen equal_or_greater = abs(ri_F) >= abs(`waldf')
	qui sum equal_or_greater, meanonly
	local waldp_ri = r(mean)
restore

** Cluster RA omnibus balance check (error adjust only)
* (Note: this is CR1, not CR2)
qui reg zcluster cov1 cov2, cluster(buildingid)
qui test cov1=cov2=0
local waldb_p = r(p)
local waldb_f = r(F)

** Repeat using randomization inference
capture program drop ri_draw
program define ri_draw, rclass
	capture drop riZ
	qui sum zcluster
	local zsum = r(sum)
	cluster_ra riZ, cluster_var(buildingid)
	qui reg riZ cov1 cov2, cluster(buildingid)
	qui test cov1=cov2=0
	return scalar ri_F = r(F)
end
preserve
	simulate ///
	ri_F = r(ri_F), ///
	reps(1000) nodots: ///
	ri_draw
	gen equal_or_greater = abs(ri_F) >= abs(`waldb_f')
	qui sum equal_or_greater, meanonly
	local waldb_p_ri = r(mean)
restore

** Cluster RA balance check: Wald test
preserve
	collapse (mean) cov1 cov2 (first) zcluster, by(buildingid)
	qui reg zcluster cov1 cov2, vce(hc2)
	qui test cov1=cov2=0
restore
local wald_clp = r(p)
local wald_clf = r(F)

** Repeat using randomization inference
capture program drop ri_draw
program define ri_draw, rclass
	capture drop riZ
	qui sum zcluster
	local zsum = r(sum)
	cluster_ra riZ, cluster_var(buildingid)
	preserve
		collapse (mean) cov1 cov2 (first) riZ, by(buildingid)
		qui reg riZ cov1 cov2, vce(hc2)
		qui test cov1=cov2=0
	restore
	return scalar ri_F = r(F)
end
preserve
	simulate ///
	ri_F = r(ri_F), ///
	reps(1000) nodots: ///
	ri_draw
	gen equal_or_greater = abs(ri_F) >= abs(`wald_clf')
	qui sum equal_or_greater, meanonly
	local wald_clp_ri = r(mean)
restore

** Blocked RA omnibus balance check
encode blockv2, gen(fblockv2)
qui reg zblockv2 cov1 cov2 i.fblockv2, vce(hc2)
test cov1=cov2=0
local wald_blp = r(p)
local wald_blf = r(F)

** Repeat using randomization inference
capture program drop ri_draw
program define ri_draw, rclass
	capture drop riZ
	block_ra riZ, block_var(blockv2)
	qui reg riZ cov1 cov2 i.fblockv2, vce(hc2)
	qui test cov1=cov2=0
	return scalar ri_F = r(F)
end
preserve
	simulate ///
	ri_F = r(ri_F), ///
	reps(1000) nodots: ///
	ri_draw
	gen equal_or_greater = abs(ri_F) >= abs(`wald_blf')
	qui sum equal_or_greater, meanonly
	local wald_blp_ri = r(mean)
restore

** Organize
matrix omnibus = J(4, 2, .)
matrix rownames omnibus = "complete RA" "clustered RA (se only)" "clustered RA" "blocked RA"
matrix colnames omnibus = "wald_p" "ri_p"
matrix omnibus[1,1] = round(`waldp', 0.001)
matrix omnibus[1,2] = round(`waldp_ri', 0.001)
matrix omnibus[2,1] = round(`waldb_p', 0.001)
matrix omnibus[2,2] = round(`waldb_p_ri', 0.001)
matrix omnibus[3,1] = round(`wald_clp', 0.001)
matrix omnibus[3,2] = round(`wald_clp_ri', 0.001)
matrix omnibus[4,1] = round(`wald_blp', 0.001)
matrix omnibus[4,2] = round(`wald_blp_ri', 0.001)

** Review
matrix list omnibus

********************** 4.7

use `restore', clear

stop

* ssc install xbalance.
* See "help xbalance" for additional Stata setup instructions.
* This is calling the RItools R package, so you will need R installed.
* The necessary path to Rterm.exe may look something like this:
global Rterm_path "C:\Program Files\R\R-4.2.1\bin\x64\Rterm.exe"

** Complete RA
gen single_block = 1 // To force only unstratified balance testing
label val zblockV2 // Remove value labels first
label val z
xbalance z single_block cov1 cov2

** Block RA
* zblockV2 instead of zblockV3
* zblockV3 is generated using blockTools, with no Stata equivalent
xbalance zblockV2 blockV2 cov1 cov2 
