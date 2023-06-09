//==============================================================================
//	STEP 4 : PROBIT REGRESSIONS (2012) 
//==============================================================================
//  Instructions:
//------------------------------------------------------------------------------
//  (1) Download all codes to <C:\ihds> folder and run from there.
//==============================================================================
//	Written	by Suresh Paul, Algorithm Basics LLC, sureshlazaruspaul@gmail.com
//	version: 05/22/2023
//==============================================================================



clear all 
capture log close 
set more off 

//------------------------------------------------------------------------------
// open logfile ... 
//------------------------------------------------------------------------------
log using "C:\ihds\step4_probit_12test", replace

//------------------------------------------------------------------------------
// set dir for data
//------------------------------------------------------------------------------
local dirpath "C:\ihds"
cd `dirpath'



//------------------------------------------------------------------------------
// Some Macros: 
//	Program append2empty: append to an empty dataset 
//	Program pe: Print Execute with loops
//------------------------------------------------------------------------------

capture program drop append2empty
program define append2empty
	tempvar qwerty // create tempvars 
		gen `qwerty' = 1
	append using `1' // append all datasets  
end 
/*------------------------------------------------------------------------------
clear
append_empty somedataset
------------------------------------------------------------------------------*/

program define pe
	if `"`0'"' != " " {
		display as text `"`0'"'
		`0'
		display("")
	}
end



/*------------------------------------------------------------------------------
 control variables used 
--------------------------------------------------------------------------------
 copc1               // ln(copc) 
 cs4x                // Distance to School (kms) 
 nchildren           // Number of children in HH
 hheduc              // Education of HH Head
 hhfemale            // Female HH head dummy 
 hhage               // Age of HH head
 mdm                 // midday meals? 
 female              // Female dummy 
 age                 // Age of child (8, 9, 10, 11) 
 brahmin forward othback dalit adivasi mcsj // caste & religion dummies
--------------------------------------------------------------------------------
 Instruments:
--------------------------------------------------------------------------------
vi4a                 // Proportion of HHs electrified in a village
tlinesdensity1       // Trasmission line length density (Line length per sqkm)
------------------------------------------------------------------------------*/



/*------------------------------------------------------------------------------
 define local macros ...
------------------------------------------------------------------------------*/
local elecvar $electricity

local depvars `" "walkfuel_fadu" "' // dependent variables ... "agematch_ed5x" "read" "math" "write"

local cvars `" $cvars "' 
local cvars1 `" $cvars1 "' 

local ovars  i.stateid i.groups8x // state, caste, religion, year dummies 

local instr `" $instrvars "' // instruments used

global sex = 0 // 0 = male, 1 = female 







timer on 1 

foreach m of local depvars { 

	pe use "$indat" , clear 

        //--------------------------------------------------------------------------
        // 1. data restrictions .. 
        //--------------------------------------------------------------------------
	pe gen sample = inlist(adolescents, 1)     & /// age between 8 and 11 
			inlist(rural      , 1)     & /// urban = 0
			inlist(ed4x       , 1)       //  currently enrolled in school

	pe tab sample, missing 
        //--------------------------------------------------------------------------


	gen wtrtym = cond(!mi(watertime_fadu), watertime_fadu, cond(!mi(watertime_f), watertime_f , cond(!mi(watertime_m), watertime_m, .)))
	gen fultym = cond(!mi(walkfuel_fadu), walkfuel_fadu, cond(!mi(walkfuel_f), walkfuel_f , cond(!mi(walkfuel_m), walkfuel_m, .)))

	egen wtrtym1 = rowmax(watertime_fadu watertime_f watertime_m) 
	egen fultym1 = rowmax(walkfuel_fadu walkfuel_f walkfuel_m)

	egen wmiss = rowmiss(watertime_fadu watertime_f watertime_m)  
	egen fmiss = rowmiss(walkfuel_fadu walkfuel_f walkfuel_m)  
	
	egen wtrtym2 = rowtotal(watertime_fadu watertime_f watertime_m) if wmiss < 3
	egen fultym2 = rowtotal(walkfuel_fadu walkfuel_f walkfuel_m) if fmiss < 3
		drop wmiss fmiss
	
	pca watertime_fadu watertime_f watertime_m
	rotate
	predict water_index
	estat kmo

        //--------------------------------------------------------------------------
        // 2. drop missing obervations of all study variables ..
        //--------------------------------------------------------------------------
	local n_instr : word count `instr'
        local inst // instruments array 
	forval i = 1 / `n_instr' {
		local inst `inst' `: word `i' of `instr''
	} 

	gen missing = 0

        local x
        local x  `m' `cvars1' `elecvar' `inst' stateid groups8x 
        foreach l of local x {
		pe replace missing = 1 if mi( `l' )
	} 

	pe replace sample = 0 if inlist(sample, 1) & inlist(missing, 1) 
		pe drop missing 
	pe tab sample, missing 
        //--------------------------------------------------------------------------





        //--------------------------------------------------------------------------
        // 3. multicollinearity check: states with 100% electrification 
        //--------------------------------------------------------------------------
        pe tabulate stateid 
        pe bys stateid: egen meanelec = mean( `elecvar' ) if inlist( sample, 1 )
        pe tabulate stateid if meanelec == 1 // which states will be dropped?

	pe replace sample = 0 if inlist(sample, 1) & inlist(meanelec, 1) 
		pe drop meanelec 
	pe tab sample, missing 
        //--------------------------------------------------------------------------






        //--------------------------------------------------------------------------
        // 4. calculate no. of PSUs in each STRATA 
        //--------------------------------------------------------------------------
	pe bys stateid distid psuid : gen uniqpsu = _n == 1 
		pe bys stateid distid : replace uniqpsu = sum(uniqpsu) 
		pe bys stateid distid : replace uniqpsu = uniqpsu[_N] 

	// must be more than one PSU in each STRATA 
	pe replace sample = 0 if inlist(sample, 1) & uniqpsu <= 1
		pe drop uniqpsu 
	pe tab sample, missing 
        //--------------------------------------------------------------------------







	pe keep if inlist( sample, 1 ) // reg sample 
		pe tab sample, missing 

        //--------------------------------------------------------------------------
        // 5. describe variables ... 
        //--------------------------------------------------------------------------
        pe mdesc `x' 
        //--------------------------------------------------------------------------




        di _n "-----------------------------------------------------------------" 
        di    "{title:PROBIT REGRESSIONS}"
        di    "-----------------------------------------------------------------" 
        di    "Education variable: `m'" 
        di    "-----------------------------------------------------------------"

        // pe probit `m' i.`elecvar' `cvars' `ovars'  , vce( robust ) 
           pe probit `m' i.`elecvar' `cvars'  , vce( robust ) 

		pe eststo: margins , dydx(`elecvar') post // marginal effects

} 


// output coefficients to Latex
esttab using step4_probit_12test.tex , b(a6) se(6) wide starlevels("\sym{*}" 0.100 "\sym{**}" 0.050 "\sym{***}" 0.010) replace 

esttab, cells(_sigsign) nogap starlevels("+/-" 0.10 "++/--" 0.05 "+++/---" 0.01)
	eststo clear


timer off 1 
timer list 1 

local min = int(r(t1)/60)  
local sec = mod(r(t1),60)

//------------------------------------------------------------------------------
// total runtime for this program : 
//------------------------------------------------------------------------------
di "total runtime for this program : `min' minutes `sec' seconds"  

timer clear 

	
// CLOSE OUTPUT
log close
