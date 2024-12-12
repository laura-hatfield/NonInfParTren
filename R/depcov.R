#' Health insurance coverage data for young adults 2008-11
#'
#' SIPP 2008 panel data on health insurance coverage of young adults originally analyzed in 
#' Antwi, Yaa Akosa, Asako S Moriya, and Kosali Simon. 
#' “Effects of Federal Policy to Insure Young Adults: Evidence from the 2010 Affordable Care Act’s Dependent-Coverage Mandate.” 
#' American Economic Journal: Economic Policy 5, no. 4 (November 2013): 1–28. https://doi.org/10.1257/pol.5.4.1.
#' This version of the data set modified for use in the NonInfParTren package by Bilinski and Hatfield.
#'
#'
#' @format A data frame with 150,996 rows and 41 variables:
#' \describe{
#'   \item{anyhi}{binary indicator of health insurance coverage from any source}
#'   \item{emphi_dep}{binary indicator of employer-sponsored health insurance coverage as a dependent}
#'   \item{emphi}{binary indicator of employer-sponsored health insurance coverage as an employee}
#'   \item{indiv}{indicator of private individual coverage}
#'   \item{govhi}{indicator of government-provided health insurance}
#'   \item{fedelig}{binary indicator of being eligible for the ACA's dependent coverage provision, i.e., age 19-25}
#'   \item{st_with_law}{binary indicator of a state with an existing dependent coverage law prior to the ACA}
#'   \item{fipstate}{FIPS code for the state}
#'   \item{trend}{integer counting up from 0 (Aug 2008) to 39 (Nov 2011); created by: (year-2008)*12 +(month-8)}
#'   \item{year}{year (2008 to 2011)}
#'   \item{month}{month of year (1 to 12)}
#'   \item{ageXX}{dummy variables for each single year of age from XX=16 to XX=29 (omitting age 26 because of ambiguity about coverage at that age)}
#'   \item{female}{binary indicator of female sex}
#'   \item{hispanic}{binary indicator of Hispanic ethnicity}
#'   \item{white}{binary indicator of White race}
#'   \item{asian}{binary indicator of Asian race}
#'   \item{other}{binary indicator of other race}
#'   \item{mar}{binary indicator of married}
#'   \item{student}{binary indicator of student}
#'   \item{fpl_ratio}{household income as a share of federal poverty line}
#'   \item{fpl_ratio_2}{household income as a share of federal poverty line, squared}
#'   \item{ue}{monthly unemployment at state level, http://127.0.0.1:21599/graphics/plot_zoom_png?width=760&height=855}
#'   \item{ue_treat}{interaction of unemployment and an indicator for treatment group}
#'   \item{weight}{person weights from SIPP, standardized to sum to 1}
#'   \item{month.id}{}
#'   \item{enact.trt.month}{an indicator for the period after ACA enactment but before implementation}
#'   \item{impl.trt.month}{an indicator for the period after ACA implementation}
#'   \item{young.ctrl}
#'   \item{old.ctrl}
#' }
#' 
#' @source https://www.openicpsr.org/openicpsr/project/114840/version/V1/view
"depcov"

# First post-period: months from enactment to just before implementation (March 2010 to Sept 2010)
# Second post-period: implementation to the end of the study (Oct 2010 to Nov 2011)
# Fedelig : age>=19 and age<26
# elig_mar10 & elig_oct10 : interactions of post-reform dummies and fedelig
# trend : (year-2008)*12 +(month-8)
effects = c("mar_sep10", "after_oct10", "fedelig", "elig_mar10", "elig_oct10")
# These are all the additional regressors;
# From the footnote of Table 2 in Akosa Antwi:
# Other regressors are 
#   an indicator for the period after ACA enactment but before implementation, 
#   an indicator for the period after ACA implementation, 
#   an indicator for each year of age, 
#   year-specific fixed effects,
#   month-specific fixed effects,
#   time trend, 
#   state fixed effects, 
#   gender, 
#   race/ethnicity, 
#   marital status, 
#   student status, 
#   household income as a share of federal poverty line and its squared term, 
#   monthly unemployment at state level, http://127.0.0.1:21599/graphics/plot_zoom_png?width=760&height=855
#   interaction of unemployment and an indicator for treatment group