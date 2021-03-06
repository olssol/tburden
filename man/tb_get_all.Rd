% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tb_inference.R
\name{tb_get_all}
\alias{tb_get_all}
\title{Overall results}
\usage{
tb_get_all(
  dat_tb,
  dat_surv,
  inx_bs = 0,
  fml_surv = "~BASE+AGE+SEX+STRATA1+P1TERTL",
  fml_tb = "~AGE+SEX+STRATA1+P1TERTL",
  imp_m = 5,
  fit_tb = TRUE,
  date_dbl = "2020-03-01",
  uti_gamma = c(Progression = 0.2, Death = 0.5),
  scenario = "scenario",
  mdl_surv = c("msm", "weibull"),
  keep_est_par = FALSE,
  ...,
  seed = NULL
)
}
\arguments{
\item{dat_tb}{tumor burden dataset}

\item{dat_surv}{survival dataset}

\item{fml_surv}{formula for survival model}

\item{fml_tb}{formula for tumor burden model}

\item{imp_m}{number of imputations for each subject}

\item{fit_tb}{whether fit tumor burden curve or only use the observed tumor
burden}

\item{date_dbl}{database lock date, i.e., analysis date}

\item{scenario}{label of scenarios}

\item{mdl_surv}{regression method for survival. msm: multi-state model;
weibull: weibull regression for PFS}

\item{uti_gamm}{utility gamma for progression and death}
}
\value{
tb_mat: tumor burden outcome for each patient
}
\description{
Overall results
}
