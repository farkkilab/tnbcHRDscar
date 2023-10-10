\name{tnbcHRDscar}
\alias{tnbcHRDscar}
\title{get.tnbcHRDscars}
\usage{
get.tnbcHRDscars(a, chrominfo ="grch38")
}
\description{
Get the number of allelic imbalances associated with HRD in triple negative breast cancer with the function get.tnbcHRDscars
}
\examples{
a <-read.table("/examples/segments.txt", header=T)
get.tnbcHRDscars(a, chrominfo ="grch38")
}
