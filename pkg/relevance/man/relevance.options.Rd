\name{relevance.options}
\alias{relevance.options}
\alias{p.symbols}
\alias{rlv.symbols}

\docType{data}
\title{Options for the relevnance Package
}
\description{List of options used in the relevnance package
  to select items and formats for printing inference elements
}
\usage{
relevance.options
rlv.symbols
p.symbols
}
\format{
  The format is:
List of 22
 $ digits.reduced        : 3
 $ testlevel             : 0.05
 $ show.confint          : TRUE
 $ termtable             : TRUE
 $ vif                   : TRUE
 $ rlvThres              : stand=0.1, rel=0.1, prop=0.1, coef=0.1, drop=0.1, pred=0.05
 $ show.ifc              : "relevance"
 $ show.ifc.relevance    : "Rle" "Rlp" "Rls" "Rls.symbol"
 $ show.ifc.test         : "Sig0" "p.symbol"
 $ show.ifc.classical    : "statistic" "p.value" "p.symbol"
 $ show.term.relevance   : "df" "R2.x" "coefRlp" "coefRls" ...
 $ show.term.test        : "df" "ciLow" "ciUp" "R2.x" ...
 $ show.term.classical   : "statistic" "df" "ciLow" "ciUp" ...
 $ show.termeff.relevance: "coef" "coefRls.symbol"
 $ show.termeff.test     : "coef" "p.symbol"
 $ show.termeff.classical: "coef" "p.symbol"
 $ show.symbollegend     : TRUE
 $ na.print              : "."
 $ p.symbols              : List, see below
 $ rlv.symbols            : List, see below

rlv.symbols
List
$ symbol  : " " "." "+" "++" "+++"
$ cutpoint: -Inf 0 1 2 5 Inf

p.symbols
List
$ symbol  : "***" "**" "*" "." " "
$ cutpoint: 0 0.001 0.01 0.05 0.1 1
}

\examples{
relevance.options
options(relevance.options) ## restores the package's default options
}
\keyword{utilities}
