#' A data set containing the clinical data of the subjects sequenced for MRSA study
#'
#' Data from a study to define the nasal microbiome of hospital inpatients who are persistently colonized with methicillin-resistant Staphylococcus aureus (MRSA) compared with matched, non-colonized controls. See original paper for matching proceedure.
#'
#' @format A 52x15 data.frame
#' \describe{
#'   \item{\code{Lib}}{A character vector of sequencing library names}
#'   \item{\code{Subject}}{Subject names}
#'   \item{\code{Age}}{Subject's Ages in years}
#'   \item{\code{Antibiotics}}{Whether or not subjects were on some form of antibiotics. 'N' = No, 'Y' = Yes}
#'   \item{\code{Aureus_Positve}}{Whether or not subjects were positive for S.Aureus (methicillin-sensitive or methicillin-resistent)}
#'   \item{\code{Case_or_Control}}{Whether or not the subject was MRSA positive (methicillin-resistent S.Aureus positive). 'Case' = MRSA positive}
#'   \item{\code{Decade}}{Subjects age in decades. '8' means the subject was in their 80s}
#'   \item{\code{Diabetes}}{Whether or not the subject had diabets. 'N' = No, 'Y' = Yes}
#'   \item{\code{Match}}{Numeric variable indicating matched subjects. See original paper for full matching procedure.}
#'   \item{\code{MRSA_Postive}}{Whether or not the subject was MRSA positive. 'Y' = MRSA positive}
#'   \item{\code{Nasal_Steroids}}{Whether or not the subject was using nasal steroids. 'N' = No, 'Y' = Yes}
#'   \item{\code{Nursing_Home}}{Whether or not the subject was staying in a nursing home. 'N' = No, 'Y' = Yes}
#'   \item{\code{Sex}}{Subject sex. All male as subjects were matched on sex}
#'   \item{\code{Smoking}}{Subject's smoking status. 'F' = Former, 'N' = Never, 'Y' = Yes/current smoker}
#'
#' }
#' @source \url{http://dx.doi.org/10.1016/j.jinf.2015.08.008}
"mrsa_clin"
