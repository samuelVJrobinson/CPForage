% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/efficiency.R
\name{efficiency}
\alias{efficiency}
\title{Efficiency}
\usage{
efficiency(L, L_max, e, d, v, h, f, l, p_i, c_i, c_f, H, alphaVal, betaVal, S)
}
\arguments{
\item{L}{Load size (\eqn{\muL})}

\item{L_max}{Maximum load size (\eqn{\muL})}

\item{e}{Energetic value of nectar (J/\eqn{\muL})}

\item{d}{Distance from hive (m)}

\item{v}{Unloaded flight speed from hive (m/s)}

\item{h}{Handling time per flower (s)}

\item{f}{Flight time between flowers (s)}

\item{l}{Maximum standing crop per flower (\eqn{\muL})}

\item{p_i}{Licking speed for nectar (\eqn{\muL/s})}

\item{c_i}{Cost of non-flying behaviour (J/s)}

\item{c_f}{Cost of flight (J/s)}

\item{H}{Time spent inside hive (s)}

\item{alphaVal}{Increase in metabolic rate with load (J/(s*\eqn{\muL}))}

\item{betaVal}{Reduction in flight speed with increase in load}

\item{S}{Competition term (0-1)}
}
\value{
Efficiency (dimensionless). \eqn{Efficiency = \frac{Gains - Foraging
 Loss - Travel Loss - Hive Loss}{Loading Loss + Travel Loss + Hive Loss}}
 Called by \code{curr_i}.
}
\description{
Returns Efficiency (ratio of profits/losses) for a forager in a given cell.
}
\examples{
efficiency(L=50,L_max=59.5,e=14.35,d=100,v=7.8,
 h=1.5,f=0.86,l=1,p_i=1,c_i=0.0042,c_f=0.05,H=100,S=0.5,
 alphaVal=0.013,betaVal=0.102/59.5)

}
