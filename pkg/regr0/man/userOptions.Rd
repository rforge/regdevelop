\name{userOptions}
\alias{userOptions}
\alias{getUserOption}
\alias{UserOptions}
\alias{UserDefault}
\title{Set and Get User "Session" Options that influence "regr0"s behavior}
\description{
  The user can set (and get) "Session" options which influence the
  behavior \pkg{regr0} functions, notably \code{\link{stamp}()}.
}
\usage{
userOptions  (x=NULL, default=NULL, list=NULL, ...)
getUserOption(x, default = NULL)
UserOptions
UserDefault
}
\arguments{
  \item{x}{option to query or set}
  \item{default}{default value}%% FIXME
  \item{list}{....}%% FIXME
  \item{...}{currently unused}
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{
FIXME
}
\author{Werner Stahel}
%% \note{
%% %%  ~~further notes~~
%% }
\seealso{
  \code{\link{stamp}}; \R's own predefined \code{\link{options}()}.
}
\examples{
## see example(stamp)
%% FIXME: add more
}
\keyword{utilities}
