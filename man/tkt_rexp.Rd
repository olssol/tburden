% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toolkits.R
\name{tkt_rexp}
\alias{tkt_rexp}
\title{Simulate time to events in days}
\usage{
tkt_rexp(
  ntot,
  hazard = NULL,
  median_mth = 5,
  annual_drop = NULL,
  mth_to_days = 30.4,
  take_floor = TRUE
)
}
\arguments{
\item{ntot}{total number of patients}

\item{hazard}{hazard for event}

\item{median_mth}{median survival in months}

\item{annual_drop}{annual drop rate}
}
\description{
Simulate time to events in days
}
