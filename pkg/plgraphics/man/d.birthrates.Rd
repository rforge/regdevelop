\name{d.birthrates}
\alias{d.birthrates}
\alias{d.birthratesVars}

\docType{data}
\title{Birthrates in Swiss Districts}
\description{Standardized fertility measure and socio-economic indicators for
  each of 182 districts of Switzerland at about 1888.
  This is an extended version of the \code{swiss} dataset of standard R.
}
\usage{
  data("d.birthrates")
  data("d.birthratesVars")
}
\format{
  \code{d.birthrates}:
  A data frame with 182 observations on the following 25 variables.
  \describe{
    \item{\code{fertility}}{Common standardizedfertility measure, see details}
    \item{\code{fertTotal}}{Alternative fertility measure}
    \item{\code{infantMort}}{Infant mortality}
    \item{\code{catholic}}{percentage of members of the catholic church}
    \item{\code{single24}}{percentage of women aged 20-24 who are single}
    \item{\code{single49}}{percentage of women aged 45-49 who are single}
    \item{\code{eAgric}}{Proportion male labor force in agriculture}
    \item{\code{eIndustry}}{Proportion male labor force in industry}
    \item{\code{eCommerce}}{Proportion male labor force in trade}
    \item{\code{eTransport}}{Proportion male labor force in transportation}
    \item{\code{eAdmin}}{Proportion male labor force in public service}
    \item{\code{german}}{percentage of German}
    \item{\code{french}}{percentage of French}
    \item{\code{italian}}{percentage of Italian}
    \item{\code{romansh}}{percentage of Romansh}
    \item{\code{gradeHigh}}{Prop. high grade in draftees exam}
    \item{\code{gradeLow}}{Propr. low grade in draftees exma}
    \item{\code{educHigh}}{Prop. draftees with > primary educ.}
    \item{\code{bornLocal}}{Proportion living in commune of birth}
    \item{\code{bornForeign}}{Proportion born in foreign country}
    \item{\code{sexratio}}{Sex ratio (M/F)}
    \item{\code{canton}}{Canton Name}
    \item{\code{district}}{District Name}
    \item{\code{altitude}}{altitude in three categories: low, medium, high}
    \item{\code{language}}{dominating language: german, french, italian,
      romansh}
  }
  \code{d.birthratesVars}:
  Data.frame that contains the descriptions of the variables just read.
}
\details{
  \code{?swiss} says:\cr
  (paraphrasing Mosteller and Tukey):\cr
  Switzerland, in 1888, was entering a period known as the
  'demographic transition'; i.e., its fertility was beginning to
  fall from the high level typical of underdeveloped countries.\cr

  The exact definition of fertility is as follows.\cr
  fertility = 100 * B_l/ sum m_i f_i, where \cr
  B_l = annual legitimate births,
  m_i = the number of married women in age interval i,
  and f_i = the fertility Hutterite women in the same age interval.\cr
  "Hutterite women" are women in a population that is known to be extremely
  fertile.\cr  
  Stillbirths are included.
}
\source{https://opr.princeton.edu/archive/pefp/switz.aspx
}
\references{see source
}
\examples{
data(d.birthrates)
## maybe str(d.birthrates) ; plot(d.birthrates) ...
}
\keyword{datasets}
