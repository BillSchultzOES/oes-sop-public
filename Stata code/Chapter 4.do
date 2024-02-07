
********************** 4.1

** Set root directory
cd "...Stata code"

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

********************** 4.4.2

/*

* Possible to do something similar in Stata, but seems like a bigger lift.
* Probably can take advantage of tools for clustered data analysis.

## Using the blockTools package
mvblocks <- block(dat1, id.vars="id", block.vars=c("cov1","cov2"), algorithm="optimal")
dat1$blocksV3 <- createBlockIDs(mvblocks, data=dat1, id.var = "id")
dat1$ZblockV3 <- labelTC(block_ra(blocks = dat1$blocksV3))

## Just show the first ten pairs
with(dat1,table(blocksV3,ZblockV3,exclude=c()))[1:10,]

--

## Using the quickblock package
distmat <- distances(dat1, dist_variables = c("cov1bin", "cov2"), id_variable = "id", normalize="mahalanobiz")
distmat[1:5,1:5]

--

quantile(as.vector(distmat), seq(0,1,.1))

--

## The caliper argument helps prevent ill-matched points
mvbigblock <- quickblock(distmat, size_constraint = 6L, caliper = 2.5)

## Look for missing points
table(mvbigblock,exclude=c()) # One point dropped due to caliper

--

dat1$blocksV4 <- mvbigblock
dat1$notblocked <- is.na(dat1$blocksV4) 
dat1$ZblockV4[dat1$notblocked==F] <- labelTC(block_ra(blocks = dat1$blocksV4))
with(dat1, table(blocksV4, ZblockV4, exclude=c()))[1:10,]

--

blockingDescEval <- dat1 %>% 
  group_by(blocksV4) %>%
  summarize(
    cov2diff = max(abs(cov2)) - min(abs(cov2)),
    cov1 = mean(cov1bin),
    count_in_block = n()
    )

blockingDescEval

*/

********************** 4.5

** Make an indicator for cluster membership
qui count
global ndat1 = r(N)
egen buildingID = seq(), from(1) to(10) block(1)
set seed 12345
cluster_ra zcluster, cluster_var(buildingID)
table zcluster buildingID

********************** 4.7

* ssc install xbalance.
* See "help xbalance" for additional setup instructions.
* This is calling the RItools R package, so you will need R installed.
* The necessary path to Rterm.exe may look something like this:
global Rterm_path "C:\Program Files\R\R-4.2.1\bin\x64\Rterm.exe"
gen single_block = 1 // To force only unstratified balance testing
label val zblockV2 // Remove value labels now
label val z

xbalance z single_block cov1 cov2

**********

xbalance zblockV2 blockV2 cov1 cov2 
