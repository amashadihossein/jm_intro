# Equi-distance sampled AR1 structure
add_cor_noise <- function(nsubj, nobs_per_subj, rho){
  
  ar1_cor <- function(n, rho) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    rho^exponent
  }
  
  MASS::mvrnorm(n = nsubj, mu = rep(0, nobs_per_subj), 
                Sigma = ar1_cor(n = nobs_per_subj, rho = rho))
}

add_baseline <- function(sim_ts, 
                         n_feature_group = 3, 
                         n_feature_per_group = 3,
                         rho_feature_group = .6){
  
  # browser()
  nsubj <- n_distinct(sim_ts$subjid)
  
  # weights/coefficient of baseline 
  gamma_coef <- rep(x = c(-1,0,1)[(1:n_feature_group) %% 3 + 1],
                    each = n_feature_per_group)
  
  sigma_block <- matrix(rho_feature_group, 
                        nrow = n_feature_per_group,
                        ncol = n_feature_per_group)
  diag(sigma_block) <- 1
  sigma <- Matrix::bdiag(replicate(3,sigma_block,simplify=FALSE))
  sigma <- Matrix::nearPD(sigma)$mat
  
  n_features <- n_feature_group * n_feature_per_group
  d <-  MASS::mvrnorm(n = nsubj, mu = rep(0, n_features), 
                      Sigma = sigma)
  gamma_x <- d %*% gamma_coef
  
  #scale
  gamma_x <- scale(gamma_x)
  gamma_x <- gamma_x[,1]/15 # 15 is chosen s.t. 3 SD is within (-.2, .2)
  
  d <- as.data.frame(d) %>%
    `colnames<-`(paste0("baseline",1:n_features)) %>%
    mutate(gamma_x = gamma_x) %>%
    mutate(subjid = unique(sim_ts$subjid)) %>%
    left_join(x = sim_ts, y = .) %>%
    relocate(subjid, all_of(paste0("baseline",1:n_features)),gamma_x)
  
  invisible(d)
}


sim_ts <- function(mu_profiles_true, 
                   nsubj_per_group = 50,
                   rho_ar1 = .6,
                   sigma_i_scale =  .1,
                   sigma_ij_scale = .05){
  
  grp_names <- unique(mu_profiles_true$grps)
  n_grp <- length(grp_names)
  nsubj <- n_grp * nsubj_per_group
  nobs_per_subj <- n_distinct(mu_profiles_true$t_obs)
  t_seq <- sort(unique(mu_profiles_true$t_obs))
  
  d <- expand.grid(t_obs = unique(mu_profiles_true$t_obs), 
                   subjid =paste0("subj",1:nsubj))
  d <- d %>% distinct(subjid) %>% 
    mutate(grps = rep(grp_names,nsubj_per_group)) %>%
    left_join(d,.)
  d <- d %>% left_join(x = ., y = mu_profiles_true) %>%
    mutate(t_obs = as.numeric(t_obs))
  
  noise <- add_cor_noise(nsubj =nsubj, nobs_per_subj = nobs_per_subj,
                         rho = rho_ar1) %>%
    data.frame() %>%
    `colnames<-`(as.character(t_seq)) %>%
    mutate(subjid = paste0("subj",1:nsubj)) %>%
    tidyr::pivot_longer(cols = - subjid, names_to = "t_obs", 
                        values_to = "noise_i") %>%
    mutate(noise_ij = rnorm(n = nrow(.))) %>%
    mutate(noise_i = sigma_i_scale * noise_i) %>%
    mutate(noise_ij = sigma_ij_scale * noise_ij) %>%
    mutate(t_obs = as.numeric(t_obs))
  
  d <- d %>% left_join(x = ., y = noise) %>%
    mutate(y_bar = mu_true + noise_i ) %>%
    mutate(y = y_bar + noise_ij) %>%
    mutate(t_obs_n = (t_obs - min(t_obs))/(max(t_obs) - min(t_obs)))
  
  invisible(d)
  
}

surv_inv <- function (t, sim_dts_i, u_i, gamma_x_i) {
  
  # browser()
  # h() is the hazard function and we assume a Weibull baseline hazard
  h <- function (s) {
    # the linear predictor from the mixed model evaluated at time s
    f <- as.vector(approxfun(x = sim_dts_i$t_obs_n,y = sim_dts_i$y_bar)(s))
    exp(log(shape_wb) + (shape_wb - 1) * log(s) +  gamma_x_i  + f * alpha )
  }
  
  # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
  log_s <- - integrate(h, lower = 0, upper = t)$value
  obj <- - log_s  + log(u_i)
  return(obj)
}