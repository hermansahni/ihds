/*-------------------------------------------------------------------------*
 |                                                                         
 |            STATA SUPPLEMENTAL SYNTAX FILE FOR ICPSR 37382
 |         INDIA HUMAN DEVELOPMENT SURVEY PANEL (IHDS, IHDS-II),
 |                            2005, 2011-2012
 |          (DATASET 0011: HOUSEHOLDS - LONG PANEL (PUBLIC-USE))
 |
 |
 | This Stata missing value recode program is provided for optional use with
 | the Stata system version of this data file as distributed by ICPSR.
 | The program replaces user-defined numeric missing values (e.g., -9)
 | with generic system missing "."  Note that Stata allows you to specify
 | up to 27 unique missing value codes.  Only variables with user-defined
 | missing values are included in this program.
 |
 | To apply the missing value recodes, users need to first open the
 | Stata data file on their system, apply the missing value recodes if
 | desired, then save a new copy of the data file with the missing values
 | applied.  Users are strongly advised to use a different filename when
 | saving the new file.
 |
 *------------------------------------------------------------------------*/

replace ID13 = . if (ID13 == 99)
replace GROUPS = . if (GROUPS == 99)


