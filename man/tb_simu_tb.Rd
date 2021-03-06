% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tb_simu_tb.R
\name{tb_simu_tb}
\alias{tb_simu_tb}
\title{Simulate Tumor Burden}
\usage{
tb_simu_tb(
  pt_surv,
  mdl_fit,
  visit_interval = 63,
  trt_effect = 0,
  tb_sig = 2,
  rand_effect = c(`(Intercept)` = 0, `poly(reg_t, 2, raw = TRUE)1` = 0, `poly(reg_t, 2,
    raw = TRUE)2` = 0),
  miss_rate = 0.5,
  ...,
  seed = NULL
)
}
\arguments{
\item{miss_rate}{missing probability at scaled time 1 (i.e., event time)}
}
\description{
Simulate Tumor Burden
}
