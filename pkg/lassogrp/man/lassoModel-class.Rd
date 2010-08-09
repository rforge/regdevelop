\name{lassoModel-class}
\docType{class}
\alias{lassoModel-class}
\alias{show,lassoModel-method}
\title{Class "lassoModel": Group Lasso Models}
\description{Objects of class "lassoModel" define link function, negative
  log-likelihood and corresponding gradient and Hessian for the model to
  be used in a Group Lasso problem.
}
\section{Objects from the Class}{
  Objects are typically created by calls of the form
  \code{\link{lassoModel}(...)}.
}
\section{Slots}{
  \describe{
    \item{\code{invlink}}{a function with arguments \code{eta}
      implementing the inverse link function.}
    \item{\code{link}}{a function with arguments \code{mu}
      implementing the link function.}
    \item{\code{nloglik}}{a function with arguments \code{y}, \code{mu} and
    \code{weights} implementing the \emph{negative}
    log-likelihood function.}
    \item{\code{ngradient}}{a function with arguments \code{x}, \code{y},
    \code{mu} and \code{weights}
    implementing the \emph{negative} gradient of the log-likelihood function.}
    \item{\code{nhessian}}{a function with arguments \code{x}, \code{mu} and
    \code{weights} implementing the \emph{negative} hessian of the
    log-likelihood function.}
    \item{\code{check}}{a function with argument \code{y} to check
      whether the response has the correct format.}
    \item{\code{name}}{a character name}
    \item{\code{comment}}{a character comment}
  }
}
\section{Methods}{
  \describe{
    \item{show}{object}
  }
}
\examples{
  getClass("lassoModel")
\dontshow{require("utils")}
  str( LogReg() )
}
\keyword{classes}
