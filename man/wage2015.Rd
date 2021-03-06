% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{wage2015}
\alias{wage2015}
\title{Wage Data}
\format{
Consists of white, non-hispanic individuals aging from 25 to 64 and
working more than 35 hours per week during at least 50 weeks of the year.
Excludes self-employed, individuals living in group quarters; individuals in
the military, agricultural or private household sectors; individuals with
inconsistent reports on earnings and employment status; individuals with
allocated or missing information in any of the variables used in the
analysis; individuals with hourly wage rate below $3. Contains 32,523
workers including 18,137 men and 14,382 women. The file contains the
following variables:
\describe{
\item{lnw}{log of hourly wages}
\item{weight}{CPS sampling weight}
\item{female}{gender indicator: 1 if female}
\item{exp1}{max(age-years of educ-7, 0)}
\item{exp2}{exp1^2/100}
\item{exp3}{exp1^3/100}
\item{exp4}{exp1^4/100}
\item{occ}{Aggregated occupation with 5 categories: managers, service,
sales, construction and production.}
\item{ind}{Aggregated industry with 12 categories: minery, construction,
manufacture, retail, transport, information, finance, professional,
education, leisure, services, public.}
\item{educ}{Education attainment with 5 categories: lhs (less than high
school graduate, years of educ < 12), hsg (high school graduate: years of
educ = 12), sc (some college: 13<=years of educ<=15), cg (college:
16<=years of educ<=17), ad (advanced degree: years of educ>=18).}
\item{ms}{Marital Status with 5 categories: married, widowed, separated,
divorced, and nevermarried.}
\item{region}{Regions with 4 categories: mw (midwest), so (south),
we (west), ne (northeast).}
}
}
\source{
U.S. March Supplement of the Current Population Survey (CPS) in 2015.
}
\usage{
wage2015
}
\description{
Wage Data
}
\keyword{datasets}
