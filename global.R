library(ggplot2)
library(ggdag)
require(visNetwork, quietly = TRUE)
library(patchwork)
library(dplyr)
library(survival)
library(survminer)
library(JM)
library(nlme)

library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("area", "patchwork")

data("pbc2")

gtsummary::theme_gtsummary_journal("jama")
theme_set(theme_bw())
custom_theme <- function(){
  theme_fn <- theme(legend.position='bottom',
                    legend.direction='horizontal')
  invisible(theme_fn)
}




