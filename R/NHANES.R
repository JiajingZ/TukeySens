#' @title NHANES
#' 
#' @description Dataset from The Third National Health and Nutrition Survey(NHANES III), 
#'     1988-1994, which contains information about socioeconomic background, medical
#'     record, dietary pattern, daily activities, and other health-related issues
#'     for respondents with high blood pressure (HBP) and of age 17 or older.
#'
#'     
#' @usage NHANES
#'      
#' @format A \code{tibble} with 2482 observations on 23 variables:
#' \describe{
#'   \item{white}{an indicator for whether the she or he is non-Hispanic white}
#'   \item{black}{an indicator for whether she or he is black}
#'   \item{hisp}{an indicator for whether she or he is Hispanic}
#'   \item{female}{an indicator for whether the respondent is female}
#'   \item{age_mo}{age (in months)}
#'   \item{hhsize}{household size}
#'   \item{edu}{number of years of education completed}
#'   \item{married}{an indicator for whether she or he is married}
#'   \item{widowed}{an indicator for whether she or he is widowed}
#'   \item{divorced}{an indicator whether she or he is divorced}
#'   \item{separated}{an indicator for whether she or he is separated}
#'   \item{income}{logged annual household income}
#'   \item{packyr}{pack years (number of packs smoked everyday multiplied by number of years smoked)}
#'   \item{bmi}{body mass index ((mass (lb)∕height (in)2) × 703)}
#'   \item{pulse}{radial pulse rate (beats/min)}
#'   \item{sodium}{sodium intake (mg)}
#'   \item{potassium}{potassium intake (mg)}
#'   \item{r_sodipota}{sodium–potassium ratio}
#'   \item{alcohol}{alcohol intake (g)}
#'   \item{insurance}{an indicator for whether she or he has health insurance}
#'   \item{together}{the frequency of meeting with friends or relatives per year}
#'   \item{ave_dbp}{average diastolic blood pressure (mmHg)}
#'   \item{trt_dbp}{an indicator for whether she or he takes two or more anti-hypertensives}
#' }
#' 
#' @source \url{https://wwwn.cdc.gov/nchs/nhanes/nhanes3/DataFiles.aspx}
#'  
#' @references{
#'  Dorie, V., M. Harada, N. B. Carnegie, and J. Hill (2016)
#'  A flexible, interpretable framework for assessing sensitivity to unmeasured confounding.
#'  \emph{Statistics in medicine 35}(20), 3453–3470
#'  }
"NHANES"