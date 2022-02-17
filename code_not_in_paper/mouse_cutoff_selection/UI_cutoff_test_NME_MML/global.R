library(data.table)
library(ggplot2)
library(plotly)
library(scales)
library(shinyjs)
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}
cutoff_dt=readRDS('../../../../downstream/output/mouse_analysis/UC_dNME_olap/regions_non_ts_only_tb_001_quant.rds')
