\name{d.birthrates}
\alias{d.birthrates}
\alias{d.birthratesVars}
\docType{data}
\title{Birthrates and socio-economic data of Switzerland, 1870
}
\description{The fertility of women (see Details) per district is
  to be modelled as a function of socio-economic variables.
}
\usage{
data("d.birthrates")
data("d.birthratesVars")}
\format{\code{d.birthrates} is
  a data frame with 182 observations on the following 25 variables.
  \describe{
\item{\code{fertTotal}}{Index of total fertility}
\item{\code{infantMort}}{Infant Mortality Rate}
\item{\code{catholic}}{percent Catholic}
\item{\code{single24}}{percent women aged 20-24 who are single}
\item{\code{single49}}{percent women aged 45-49 who are single}
\item{\code{eAgric}}{Proportion male labor force in agriculture}
\item{\code{eIndustry}}{Proportion male labor force in industry}
\item{\code{eCommerce}}{Proportion male labor force in trade}
\item{\code{eTransport}}{Proportion male labor force in transportation}
\item{\code{eAdmin}}{Proportion male labor force in public service}
\item{\code{german}}{percent German}
\item{\code{french}}{percent French}
\item{\code{italian}}{percent Italian}
\item{\code{romansh}}{percent Romansh}
\item{\code{gradeHigh}}{Prop. high grade in draftees exam}
\item{\code{gradeLow}}{Prop. low grade in draftees exma}
\item{\code{educHigh}}{Prop. draftees with > primary educ.}
\item{\code{bornLocal}}{Proportion living in commune of birth}
\item{\code{bornForeign}}{Proportion born in foreign country}
\item{\code{sexratio}}{Sex ratio (M/F)}
\item{\code{canton}}{Canton Name}
%%    \item{\code{canton}}{a factor with levels \code{App~Ar}
%%      \code{App~Ir} \code{Argovie} \code{B.~Stadt} \code{Basellnd}
%%      \code{Berne} \code{Fribourg} \code{Geneve} \code{Grisons}
%%      \code{Lucerne} \code{Neuchatl} \code{Schfhaus} \code{Schwyz}
%%      \code{Solothrn} \code{St~Gall} \code{Thurgove} \code{Ticino}
%%      \code{Uri} \code{Valais} \code{Vaud} \code{Zurich}} 
\item{\code{district}}{District Name}
\item{\code{altitude}}{altitude in three categories: low, medium, high}
\item{\code{language}}{dominating language: german, french, italian, romansh}
  }
    \code{d.birthratesVars} contains these descriptions.
}
\details{Exact definition of fertility:
fertility = 100 * B_l/ sum m_i f_i
where B_l = annual legitimate births, m_i = the number of married women
in age interval i, and f_i = the fertility Hutterite women
in the same age interval. Stillbirths are included.
"Hutterite women" are women in a population that is known to be extremely
fertile. 
}
\source{https://opr.princeton.edu/archive/pefp/switz.aspx
}
\references{
Van der Walle, Francine (1980).
"Education and the Demographic Transition in Switzerland."
Population and Development Review, vol. 6, no. 3.
}
\examples{
data(d.birthrates)
plot(fertility ~ eAgric, data=d.birthrates)
r.lm <- lm(fertility ~ catholic + single24 + single49 +
  eAgric + eIndustry + eCommerce + eTransport + eAdmin +
  german + french + italian + romansh + gradeHigh + gradeLow +
  educHigh + bornLocal + bornForeign + sexratio, data=d.birthrates)
r.st <- step(r.lm, k=4, trace=FALSE)
summary(r.st)
}
\keyword{datasets}
