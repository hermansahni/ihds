//==============================================================================
//	STEP 6 : BIVARIATE PROBIT REGRESSIONS (2012) 
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
log using "C:\ihds\step6_biprobit_12test", replace

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

local depvars `" "agematch_ed5x" "read" "math" "write" "' // dependent variables ... 

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

	gen elecqual3h = elecqual3 == 2
		replace elecqual3h = . if missing( elecqual3 )	
	gen elecqual3l = elecqual3 == 1
		replace elecqual3l = . if missing( elecqual3 )	
	
	gen elecqual4h = elecqual4 == 2
		replace elecqual4h = . if missing( elecqual4 )	
	gen elecqual4l = elecqual4 == 1
		replace elecqual4l = . if missing( elecqual4 )	
		
	gen elecqual5h = elecqual5 == 2
		replace elecqual5h = . if missing( elecqual5 )	
	gen elecqual5l = elecqual5 == 1
		replace elecqual5l = . if missing( elecqual5 )	


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
		di    "{title:BIVARIATE PROBIT REGRESSIONS}"
        di    "-----------------------------------------------------------------" 
        di    "Education variable: `m'" 
        di    "-----------------------------------------------------------------"

	pe cap nois biprobit (`m' = i.`elecvar' `cvars' `ovars' ) ( `elecvar' = `cvars' `ovars' `inst') , vce(robust) 

		if _rc continue

	// AME Margins (check) 
	di _n "{p}margins, at(`elecvar' = (0 1)) predict(pmarg1) force contrast(atcontrast(rb1._at)).{p_end}" 
	margins, at(`elecvar' = (0 1)) predict(pmarg1) force contrast(atcontrast(rb1._at)) 

	// AME by Hand # 1 
	di _n "{p}Margins: AME by Hand #1.{p_end}"
	predict double xb2, xb2 // linear index for equation 2
	ren `elecvar' T`elecvar'

	gen `elecvar' = 0

	predict double p0, pmarg1  // success for equation 1 with `elecvar' == 0
	predict double xb0, xb1    // index for equation 1 with `elecvar' == 0
	replace `elecvar' = 1

	predict double p1, pmarg1 // success for equation 1 with `elecvar' == 1
	predict double xb1, xb1   // index for equation 1 with `elecvar' == 1
	gen double dp=p1-p0       // calculate diff in
	sum dp

	// AME by Hand # 2
	di _n "{p}Margins: AME by Hand #2.{p_end}"
	gen double pdx=(binormal(xb1,xb2,e(rho))-binormal(xb0,xb2,e(rho)))/normal(xb2) if T`elecvar' == 1
	qui replace pdx=normal(xb1)-normal(xb0)
	su pdx

	// AME Margins (all variables) 
	di _n "{p}margins, dydx(`elecvar') predict(pmarg1) force.{p_end}"
	eststo: margins, dydx(`elecvar') predict(pmarg1) force post 

} 

	
// output coefficients to Latex
esttab using step6_biprobit_12test.tex , b(a6) se(6) wide starlevels("\sym{*}" 0.100 "\sym{**}" 0.050 "\sym{***}" 0.010) replace 

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
