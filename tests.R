library(testthat)
library(MagmaClustR)

test_check("MagmaClustR")


############ test gradients  of elbos and likelihoods #########




test_that("gradient of gr_GP() works for Squared Exponential kernel", {
  db <- tibble::tibble(Input = 1:5, Output = 2:6, Covariate = 3:7,
                       Reference = paste(Input, Covariate, sep = ':'))
  mean <- rep(0, 5)
  hp <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5)
  new_cov <- kern_to_cov(db %>% dplyr::select(- Output), "SE", hp)

  hp_v <- tibble::tibble(se_variance = 1 + 10^(-8), se_lengthscale = 0.5)
  hp_l <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5 + 10^(-8))

  deriv_v <- gr_GP(hp, db, mean, "SE", new_cov, 0)[["se_variance"]]
  deriv_l <- gr_GP(hp, db, mean, "SE", new_cov, 0)[["se_lengthscale"]]

  emp_deriv_v <- (logL_GP(hp_v, db, mean, "SE", new_cov, 0) -
                    logL_GP(hp, db, mean, "SE", new_cov, 0)) / 10^(-8)
  emp_deriv_l <- (logL_GP(hp_l, db, mean, "SE", new_cov, 0) -
                    logL_GP(hp, db, mean, "SE", new_cov, 0)) / 10^(-8)

  round(deriv_v, 3) %>% expect_equal(round(emp_deriv_v, 3))
  round(deriv_l, 3) %>% expect_equal(round(emp_deriv_l, 3))
})

## TODO: test for gr_GP() for the RQ and PERIOD kernels

test_that("gradient of logL_GP_mod() works for Squared Exponential kernel", {
  db <- tibble::tibble(Input = 1:5, Output = 2:6, Covariate = 3:7,
                       Reference = paste(Input, Covariate, sep = ':'))

  mean <- rep(0, 5)
  hp <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5)
  new_cov <- kern_to_cov(db %>% dplyr::select(- Output), "SE", hp)

  hp_v <- tibble::tibble(se_variance = 1 + 10^(-8), se_lengthscale = 0.5)
  hp_l <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5 + 10^(-8))

  deriv_v <- gr_GP_mod(hp, db, mean, "SE", new_cov, 0)[["se_variance"]]
  deriv_l <- gr_GP_mod(hp, db, mean, "SE", new_cov, 0)[["se_lengthscale"]]

  emp_deriv_v <- (logL_GP_mod(hp_v, db, mean, "SE", new_cov, 0) -
                    logL_GP_mod(hp, db, mean, "SE", new_cov, 0)) / 10^(-8)
  emp_deriv_l <- (logL_GP_mod(hp_l, db, mean, "SE", new_cov, 0) -
                    logL_GP_mod(hp, db, mean, "SE", new_cov, 0)) / 10^(-8)

  round(deriv_v, 3) %>% expect_equal(round(emp_deriv_v, 3))
  round(deriv_l, 3) %>% expect_equal(round(emp_deriv_l, 3))
})

## TODO: test for gr_GP_mod() for the RQ and PERIOD kernels

test_that("gradient of logL_GP_mod_common_hp() works", {
  db <- tibble::tibble(
    ID = rep(1:5, each = 4),
    Output = 1:20,
    Input = 2:21,
    Covariate = c(1:10, 23, 77, 1:8),
    Reference = paste(Input, Covariate, sep = ':')
  )
  mean <- tibble::tibble("Input" = db$Input, "Covariate" = db$Covariate,
                         "Reference" = db$Reference, "Output" = 0)
  hp <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5)
  new_cov <- kern_to_cov(db %>% dplyr::select(- Output), "SE", hp)

  hp_v <- tibble::tibble(se_variance = 1 + 10^(-8), se_lengthscale = 0.5)
  hp_l <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5 + 10^(-8))

  deriv_v <- gr_GP_mod_common_hp(
    hp, db, mean,
    "SE", new_cov, 0.1
  )[["se_variance"]]
  deriv_l <- gr_GP_mod_common_hp(
    hp, db, mean,
    "SE", new_cov, 0
  )[["se_lengthscale"]]

  emp_deriv_v <- (logL_GP_mod_common_hp(hp_v, db, mean, "SE", new_cov, 0.1) -
                    logL_GP_mod_common_hp(hp, db, mean, "SE", new_cov, 0.1)) / 10^(-8)

  emp_deriv_l <- (logL_GP_mod_common_hp(hp_l, db, mean, "SE", new_cov, 0) -
                    logL_GP_mod_common_hp(hp, db, mean, "SE", new_cov, 0)) / 10^(-8)

  round(deriv_v, 3) %>% expect_equal(round(emp_deriv_v, 3))
  round(deriv_l, 3) %>% expect_equal(round(emp_deriv_l, 3))
})

## TODO: test for gr_GP_mod_common_hp() for the RQ and PERIOD kernels

test_that("gradient of gr_sum_logL_GP_clust() works for SE kernel", {
  db <- tibble::tibble(Input = 1:5, Output = 2:6,
                       Covariate = 3:7, Reference = paste(1:5, 3:7, sep = ':'))
  mean <- list(
    "K1" = tibble::tibble("Input" = db$Input, "Covariate" = db$Covariate,
                          "Reference" = db$Reference, "Output" = 0),
    "K2" = tibble::tibble("Input" = db$Input, "Covariate" = db$Covariate,
                          "Reference" = db$Reference, "Output" = 0)
  )
  hp <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5)
  new_cov <- list(
    "K1" = kern_to_cov(db %>% dplyr::select(- Output), "SE", hp),
    "K2" = kern_to_cov(db %>% dplyr::select(- Output), "SE", hp)
  )
  mixture <- tibble::tibble("K1" = 0.4, "K2" = 0.6)
  hp_v <- tibble::tibble(se_variance = 1 + 10^(-8), se_lengthscale = 0.5)
  hp_l <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5 + 10^(-8))

  deriv_v <- gr_sum_logL_GP_clust(
    hp, db, mixture, mean,
    "SE", new_cov, 0
  )[["se_variance"]]
  deriv_l <- gr_sum_logL_GP_clust(
    hp, db, mixture, mean,
    "SE", new_cov, 0
  )[["se_lengthscale"]]

  emp_deriv_v <- (sum_logL_GP_clust(
    hp_v, db, mixture, mean,
    "SE", new_cov, NULL, 0
  ) -
    sum_logL_GP_clust(hp, db, mixture, mean, "SE", new_cov, NULL, 0)) / 10^(-8)
  emp_deriv_l <- (sum_logL_GP_clust(
    hp_l, db, mixture, mean,
    "SE", new_cov, NULL, 0
  ) -
    sum_logL_GP_clust(hp, db, mixture, mean, "SE", new_cov, NULL, 0)) / 10^(-8)

  round(deriv_v, 3) %>% expect_equal(round(emp_deriv_v, 3))
  round(deriv_l, 3) %>% expect_equal(round(emp_deriv_l, 3))
})


test_that("gradient of logL_GP_mod() works for Squared Exponential kernel", {
  #db <- tibble::tibble(Input = 1:5, Output = 2:6, Covariate = 3:7,
  #                                                Reference = paste(Input, Covariate, sep = ':'))

  #mean <- rep(0, 5)
  db <- simon
  hp <- tibble::tibble(ID = simon$ID,se_variance = rep(1,length(simon$ID)),
                       se_lengthscale = rep(0.5,length(simon$ID)))
  #new_cov <- kern_to_cov(db %>% dplyr::select(- Output), "SE", hp)
  new_cov <- simone$hyperpost$cov
  hp_v <- tibble::tibble(ID = simon$ID,se_variance = rep(1 + 10^(-8),length(simon$ID)),
                         se_lengthscale = rep(0.5,length(simon$ID)))
  hp_l <- tibble::tibble(ID = simon$ID,se_variance = rep(1,length(simon$ID)),
                         se_lengthscale = rep(0.5 + 10^(-8),length(simon$ID)))

  deriv_v <- gr_clust_multi_GP_common_hp_i(hp, db,new_cov, "SE", 0)[["se_variance"]]
  deriv_l <- gr_clust_multi_GP_common_hp_i(hp, db, new_cov, "SE", 0)[["se_lengthscale"]]

  emp_deriv_v <- (elbo_clust_multi_GP_common_hp_i(hp_v, db, new_cov, "SE",  0) -
                    +                         elbo_clust_multi_GP_common_hp_i(hp, db, new_cov, "SE", 0)) / 10^(-8)
  emp_deriv_l <- (elbo_clust_multi_GP_common_hp_i(hp_l, db, new_cov,"SE", 0) -
                    +                         elbo_clust_multi_GP_common_hp_i(hp, db, new_cov, "SE",0)) / 10^(-8)

  round(deriv_v, 3) %>% expect_equal(round(emp_deriv_v, 3))
  round(deriv_l, 3) %>% expect_equal(round(emp_deriv_l, 3))
})
##############################################################"
test_that("gradient of logL_GP_mod() works for Squared Exponential kernel", {
  db <- tibble::tibble(Input = 1:5, Output = 2:6, Covariate = 3:7,
                       Reference = paste(Input, Covariate, sep = ':'))

  names_col <- simon %>%
    dplyr::select(-c(.data$ID,.data$Output)) %>%
    names()

  ## Keep 6 significant digits for entries to avoid numerical errors and
  ## Add a Reference column for identification and sort according to it
  simon <- simon %>% purrr::modify_at(tidyselect::all_of(names_col),signif) %>%
    tidyr::unite("Reference",
                 tidyselect::all_of(names_col),
                 sep=":",
                 remove = FALSE) %>%
    tidyr::drop_na() %>%
    dplyr::group_by(.data$ID) %>%
    dplyr::arrange(.data$Reference, .by_group = TRUE) %>%
    dplyr::ungroup()
  all_input <- simon %>%
    dplyr::pull(.data$Reference) %>%
    unique() %>%
    sort()
  m_k <- list(); ID_k <- paste0("K", 1:3)
  ## Create a list named by cluster with evaluation of the mean at all Input
  for (k in 1:3) {
    m_k[[ID_k[k]]] <- rep(0, length(all_input))}
  ## Create a list named by cluster with evaluation of the mean at all Input
  for (k in 1:3) {
    m_k[[ID_k[k]]] <- rep(0, 20)}
  hp <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5)
  new_cov <- kern_to_cov(db %>% dplyr::select(- Output), "SE", hp)

  hp_v <- tibble::tibble(se_variance = 1 + 10^(-8), se_lengthscale = 0.5)
  hp_l <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5 + 10^(-8))

  deriv_v <- gr_GP_mod_common_hp_k(hp, db= simone$hyperpost$mu, m_k, "SE", post_cov = NULL, 1e-10,categorial = TRUE)[["se_variance"]]
  deriv_l <- gr_GP_mod_common_hp_k(hp, db= simone$hyperpost$mu, m_k, "SE", post_cov = NULL, 1e-10,categorial = TRUE)[["se_lengthscale"]]

  emp_deriv_v <- (elbo_GP_mod_common_hp_k(hp_v, db= simone$hyperpost$mu, m_k, "SE", post_cov = NULL, 1e-10,categorial = TRUE) -
                    elbo_GP_mod_common_hp_k(hp, db= simone$hyperpost$mu, m_k, "SE", post_cov = NULL, 1e-10,categorial = TRUE)) / 10^(-8)
  emp_deriv_l <- (elbo_GP_mod_common_hp_k(hp_l, db= simone$hyperpost$mu, m_k, "SE", post_cov = NULL, 1e-10,categorial = TRUE) -
                    elbo_GP_mod_common_hp_k(hp, db= simone$hyperpost$mu, m_k, "SE", post_cov = NULL, 1e-10,categorial = TRUE)) / 10^(-8)

  round(deriv_v, 3) %>% expect_equal(round(emp_deriv_v, 3))
  round(deriv_l, 3) %>% expect_equal(round(emp_deriv_l, 3))
})
##############################################################"
test_that("gradient of elbo_clust_multi_GP_common_hp_i() works for Squared Exponential kernel", {
  db <- tibble::tibble(Input = 1:5, Output = 2:6, Covariate = 3:7,
                       Reference = paste(Input, Covariate, sep = ':'))

  names_col <- simon %>%
    dplyr::select(-c(.data$ID,.data$Output)) %>%
    names()

  ## Keep 6 significant digits for entries to avoid numerical errors and
  ## Add a Reference column for identification and sort according to it
  simon <- simon %>% purrr::modify_at(tidyselect::all_of(names_col),signif) %>%
    tidyr::unite("Reference",
                 tidyselect::all_of(names_col),
                 sep=":",
                 remove = FALSE) %>%
    tidyr::drop_na() %>%
    dplyr::group_by(.data$ID) %>%
    dplyr::arrange(.data$Reference, .by_group = TRUE) %>%
    dplyr::ungroup()
  all_input <- simon %>%
    dplyr::pull(.data$Reference) %>%
    unique() %>%
    sort()

  hp <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5)
  new_cov <- kern_to_cov(db %>% dplyr::select(- Output), "SE", hp)

  hp_v <- tibble::tibble(se_variance = 1 + 10^(-6), se_lengthscale = 0.5)
  hp_l <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5 + 10^(-6))

  deriv_v <- gr_clust_multi_GP_common_hp_i(hp,db=data,hyperpost = NULL,"SE",1e-10,latents = simone$hyperpost[c("Z","mu","y_star")],categorial = TRUE)[["se_variance"]]
  deriv_l <- gr_clust_multi_GP_common_hp_i(hp,db=data,hyperpost = NULL,"SE",1e-10,latents = simone$hyperpost[c("Z","mu","y_star")],categorial = TRUE)[["se_lengthscale"]]

  emp_deriv_v <- ( elbo_clust_multi_GP_common_hp_i(hp_v,db=data,hyperpost = NULL,"SE",1e-10,latents = simone$hyperpost[c("Z","mu","y_star")],categorial = TRUE) -
                     elbo_clust_multi_GP_common_hp_i(hp,db=data,hyperpost = NULL,"SE",1e-10,latents = simone$hyperpost[c("Z","mu","y_star")],categorial = TRUE)) / 10^(-6)
  emp_deriv_l <- ( elbo_clust_multi_GP_common_hp_i(hp_l,db=data,hyperpost = NULL,"SE",1e-10,latents = simone$hyperpost[c("Z","mu","y_star")],categorial = TRUE) -
                     elbo_clust_multi_GP_common_hp_i(hp,db=data,hyperpost = NULL,"SE",1e-10,latents = simone$hyperpost[c("Z","mu","y_star")],categorial = TRUE)) / 10^(-6)

  round(deriv_v, 3) %>% expect_equal(round(emp_deriv_v, 3))
  round(deriv_l, 3) %>% expect_equal(round(emp_deriv_l, 3))
})
######################### test kern_to_cov_inv #####################

test_that("kern_to_cov() works for scalar inputs", {
  hp <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  input <- c(2, 3, 4)

  res <- matrix(NA, ncol = 3, nrow = 3)
  for (i in 1:3)
  {
    for (j in 1:3) {
      res[i, j] <- se_kernel(input[i], input[j], hp)
    }
  }
  res <- res %>%
    `rownames<-`(as.character(input)) %>%
    `colnames<-`(as.character(input))

  kern_to_cov(input, "SE", hp) %>% expect_equal(res)
})

test_that("kern_to_cov() works for vector inputs", {
  hp <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  input <- data.frame(Input = c(1, 2, 3), Cov1 = c(2, 3, 4))

  res <- matrix(NA, ncol = 3, nrow = 3)
  for (i in 1:3)
  {
    for (j in 1:3) {
      res[i, j] <- se_kernel(input[i, ], input[j, ], hp)
    }
  }

  ref = paste(input$Input, Cov1 = input$Cov1, sep=":")
  res <- res %>%
    `rownames<-`(ref) %>%
    `colnames<-`(ref)

  kern_to_cov(input, "SE", hp) %>% expect_equal(res)
})

test_that("matrix, dataframe and tibble work the same", {
  hp <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  df <- data.frame(Input = c(1, 2, 3), Cov1 = c(2, 3, 4))
  tib <- df %>% tibble::as_tibble()
  mat <- df %>% as.matrix()

  kern_to_cov(df, "SE", hp) %>% expect_equal(kern_to_cov(mat, "SE", hp))
  kern_to_cov(df, "SE", hp) %>% expect_equal(kern_to_cov(tib, "SE", hp))
  kern_to_cov(tib, "SE", hp) %>% expect_equal(kern_to_cov(mat, "SE", hp))
})

test_that("1D-matrix and vector work the same", {
  hp <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  vec <- c(1, 2, 3)
  mat <- vec %>% as.matrix()

  kern_to_cov(vec, "SE", hp) %>% expect_equal(kern_to_cov(mat, "SE", hp))
})

test_that("dimension names are correct", {
  hp <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  df <- data.frame(Input = c(5, 6, 7), Cov1 = c(2, 3, 4),
                   Reference = c('5:2', '6:3', '7:4'))
  df2 <- data.frame(Cov1 = c(2, 3, 4), Reference = c('5:2', '6:3', '7:4'),
                    Input = c(5, 6, 7))
  df3 <- data.frame(c(5, 6, 7), c(2, 3, 4))
  df4 <- data.frame(fu = c(5, 6, 7), blob = c(2, 3, 4))

  dimnames(kern_to_cov(df, "SE", hp)) %>%
    expect_equal(dimnames(kern_to_cov(df2, "SE", hp)))
  dimnames(kern_to_cov(df, "SE", hp)) %>%
    expect_equal(dimnames(kern_to_cov(df3, "SE", hp)))
  dimnames(kern_to_cov(df, "SE", hp)) %>%
    expect_equal(dimnames(kern_to_cov(df4, "SE", hp)))
})

test_that("kern_to_cov() works for custom kernels", {
  hp_se <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  hp_perio <- tibble::tibble(
    perio_variance = 2, perio_lengthscale = 1,
    period = 1
  )
  hp_rq <- tibble::tibble(rq_variance = 2, rq_lengthscale = 1, rq_scale = 1)
  df <- data.frame(Input = c(5, 6, 7), Cov1 = c(2, 3, 4))

  kern_to_cov(df, "SE", hp_se) %>%
    expect_equal(kern_to_cov(df, se_kernel, hp_se))
  kern_to_cov(df, "RQ", hp_rq) %>%
    expect_equal(kern_to_cov(df, rq_kernel, hp_rq))
  kern_to_cov(df, "PERIO", hp_perio) %>%
    expect_equal(kern_to_cov(df, perio_kernel, hp_perio))
})

test_that("kern_to_cov() works for derivative matrices", {
  hp_se <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  hp_perio <- tibble::tibble(
    perio_variance = 2, perio_lengthscale = 1,
    period = 1
  )
  hp_rq <- tibble::tibble(rq_variance = 2, rq_lengthscale = 1, rq_scale = 1)
  df <- data.frame(Input = c(5, 6, 7), Cov1 = c(2, 3, 4))

  kern_to_cov(df, "SE", hp_se) %>%
    expect_equal(kern_to_cov(df, "SE", hp_se, "se_variance"))
  kern_to_cov(df, "RQ", hp_rq) %>%
    expect_equal(kern_to_cov(df, "RQ", hp_rq, "rq_variance"))
  kern_to_cov(df, "PERIO", hp_perio) %>%
    expect_equal(kern_to_cov(df, "PERIO", hp_perio, "perio_variance"))
  ## Test for custom kernel
  kern_to_cov(df, "PERIO", hp_perio) %>%
    expect_equal(kern_to_cov(df, perio_kernel, hp_perio, "perio_variance"))
})

test_that("kern_to_cov() works for compound kernels", {
  hp_se <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  hp_perio <- tibble::tibble(
    perio_variance = 2, perio_lengthscale = 1,
    period = 1
  )
  hp_rq <- tibble::tibble(rq_variance = 2, rq_lengthscale = 1, rq_scale = 1)
  hp <- hp_se %>%
    dplyr::bind_cols(hp_perio) %>%
    dplyr::bind_cols(hp_rq)
  df <- data.frame(Input = c(5, 6, 7), Cov1 = c(2, 3, 4))

  kern_to_cov(df, "SE + RQ + PERIO", hp) %>%
    expect_equal(kern_to_cov(df, "SE", hp_se) +
                   kern_to_cov(df, "RQ", hp_rq) +
                   kern_to_cov(df, "PERIO", hp_perio))

  kern_to_cov(df, "SE * RQ * PERIO", hp) %>%
    expect_equal(kern_to_cov(df, "SE", hp_se) *
                   kern_to_cov(df, "RQ", hp_rq) *
                   kern_to_cov(df, "PERIO", hp_perio))

  kern_to_cov(df, "SE * RQ + PERIO", hp) %>%
    expect_equal(kern_to_cov(df, "SE", hp_se) *
                   kern_to_cov(df, "RQ", hp_rq) +
                   kern_to_cov(df, "PERIO", hp_perio))
})

test_that("kern_to_cov() works for compound kernels' derivatives", {
  hp_se <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  hp_perio <- tibble::tibble(
    perio_variance = 2, perio_lengthscale = 1,
    period = 1
  )
  hp_rq <- tibble::tibble(rq_variance = 2, rq_lengthscale = 1, rq_scale = 1)
  hp_lin <- tibble::tibble(lin_slope = 2, lin_offset = 1)
  hp <- hp_se %>% dplyr::bind_cols(hp_perio, hp_rq, hp_lin)

  df <- data.frame(Input = c(5, 6, 7), Cov1 = c(2, 3, 4))

  kern_to_cov(df, "SE + RQ + PERIO", hp, "se_variance") %>%
    expect_equal(kern_to_cov(df, "SE", hp_se, "se_variance"))

  kern_to_cov(df, "SE * RQ * PERIO * LIN", hp, "rq_lengthscale") %>%
    expect_equal(kern_to_cov(df, "SE", hp_se, ) *
                   kern_to_cov(df, "RQ", hp_rq, "rq_lengthscale") *
                   kern_to_cov(df, "PERIO", hp_perio) *
                   kern_to_cov(df, "LIN", hp_lin))

  kern_to_cov(df, "SE * RQ + PERIO", hp, "period") %>%
    expect_equal(kern_to_cov(df, "PERIO", hp_perio, "period"))

  kern_to_cov(df, "LIN", hp, "lin_offset") %>%
    expect_equal(kern_to_cov(df, "LIN", hp, NULL) -
                   kern_to_cov(df, "LIN", hp, "lin_slope"))
})
############### test kern_to_inv ##########
test_that("kern_to_inv() works for scalar inputs", {
  hp <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  input <- c(2, 3, 4)

  res <- matrix(NA, ncol = 3, nrow = 3)
  for (i in 1:3)
  {
    for (j in 1:3) {
      res[i, j] <- se_kernel(input[i], input[j], hp)
    }
  }
  res <- res %>%
    solve() %>%
    `rownames<-`(as.character(input)) %>%
    `colnames<-`(as.character(input))
  expect_equal(kern_to_inv(input, "SE", hp), res)
})

test_that("kern_to_inv() works for vector inputs", {
  hp <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  input <- data.frame(Input = c(1, 2, 3), Cov1 = c(2, 3, 4))

  res <- matrix(NA, ncol = 3, nrow = 3)
  for (i in 1:3)
  {
    for (j in 1:3) {
      res[i, j] <- se_kernel(input[i, ], input[j, ], hp)
    }
  }
  ref = paste(input$Input, Cov1 = input$Cov1, sep=":")

  res <- res %>%
    solve() %>%
    `rownames<-`(ref) %>%
    `colnames<-`(ref)

  kern_to_inv(input, "SE", hp) %>% expect_equal(res)
})

test_that("dimension names are correct", {
  hp <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  df <- data.frame(Input = c(5, 6, 7), Cov1 = c(2, 3, 4),
                   Reference = c('5:2', '6:3', '7:4'))
  df2 <- data.frame(Cov1 = c(2, 3, 4), Reference = c('5:2', '6:3', '7:4'),
                    Input = c(5, 6, 7))
  df3 <- data.frame(c(5, 6, 7), c(2, 3, 4))
  df4 <- data.frame(fu = c(5, 6, 7), blob = c(2, 3, 4))

  dimnames(kern_to_inv(df, "SE", hp)) %>%
    expect_equal(list(c('5:2', '6:3', '7:4'), c('5:2', '6:3', '7:4')))
  dimnames(kern_to_inv(df, "SE", hp)) %>%
    expect_equal(dimnames(kern_to_inv(df2, "SE", hp)))
  dimnames(kern_to_inv(df, "SE", hp)) %>%
    expect_equal(dimnames(kern_to_inv(df3, "SE", hp)))
  dimnames(kern_to_inv(df, "SE", hp)) %>%
    expect_equal(dimnames(kern_to_inv(df4, "SE", hp)))
})
########## test kernels #########
test_that("Squared Exponential kernel works for null distance", {
  hp <- tibble::tibble(se_variance = 1, se_lengthscale = 1)

  se_kernel(2, 2, hp)[1] %>% expect_equal(exp(1))
})

test_that("Periodic kernel works for null distance", {
  hp <- tibble::tibble(perio_variance = 1, perio_lengthscale = 1, period = 1)

  perio_kernel(2, 2, hp)[1] %>% expect_equal(exp(1))
})

test_that("Rational quadratic kernel works for null distance", {
  hp <- tibble::tibble(rq_variance = 1, rq_lengthscale = 1, rq_scale = 1)

  rq_kernel(2, 2, hp)[1] %>% expect_equal(exp(1))
})

test_that("Linear kernel works for null distance", {
  hp <- tibble::tibble(lin_slope = 1, lin_offset = 1)

  lin_kernel(0, 2, hp)[1] %>% expect_equal(exp(1))
})

test_that("gradients for the Squared Exponential kernel are valid", {
  hp <- tibble::tibble(se_variance = 1, se_lengthscale = 1)
  hp_v <- tibble::tibble(se_variance = 1 + 10^(-8), se_lengthscale = 1)
  hp_l <- tibble::tibble(se_variance = 1, se_lengthscale = 1 + 10^(-8))

  deriv_v <- se_kernel(c(1, 2), c(2, 3), hp, "se_variance")
  deriv_l <- se_kernel(c(1, 2), c(2, 3), hp, "se_lengthscale")

  emp_deriv_v <- (se_kernel(c(1, 2), c(2, 3), hp_v)[1] -
                    se_kernel(c(1, 2), c(2, 3), hp)[1]) / 10^(-8)
  emp_deriv_l <- (se_kernel(c(1, 2), c(2, 3), hp_l)[1] -
                    se_kernel(c(1, 2), c(2, 3), hp)[1]) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
})

test_that("gradients for the Periodic kernel are valid", {
  hp <- tibble::tibble(
    perio_variance = 1,
    perio_lengthscale = 1, period = pi
  )
  hp_v <- tibble::tibble(
    perio_variance = 1 + 10^(-8),
    perio_lengthscale = 1, period = pi
  )
  hp_l <- tibble::tibble(
    perio_variance = 1,
    perio_lengthscale = 1 + 10^(-8), period = pi
  )
  hp_p <- tibble::tibble(
    perio_variance = 1,
    perio_lengthscale = 1, period = pi + 10^(-8)
  )

  deriv_v <- perio_kernel(c(1, 2), c(2, 3), hp, "perio_variance")
  deriv_l <- perio_kernel(c(1, 2), c(2, 3), hp, "perio_lengthscale")
  deriv_p <- perio_kernel(c(1, 2), c(2, 3), hp, "period")

  emp_deriv_v <- (perio_kernel(c(1, 2), c(2, 3), hp_v)[1] -
                    perio_kernel(c(1, 2), c(2, 3), hp)[1]) / 10^(-8)
  emp_deriv_l <- (perio_kernel(c(1, 2), c(2, 3), hp_l)[1] -
                    perio_kernel(c(1, 2), c(2, 3), hp)[1]) / 10^(-8)
  emp_deriv_p <- (perio_kernel(c(1, 2), c(2, 3), hp_p)[1] -
                    perio_kernel(c(1, 2), c(2, 3), hp)[1]) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
  round(deriv_p, 4) %>% expect_equal(round(emp_deriv_p, 4))
})

test_that("gradients for the Rational Quadratic kernel are valid", {
  hp <- tibble::tibble(rq_variance = 1, rq_lengthscale = 1, rq_scale = 1)
  hp_v <- tibble::tibble(
    rq_variance = 1 + 10^(-8),
    rq_lengthscale = 1, rq_scale = 1
  )
  hp_l <- tibble::tibble(
    rq_variance = 1,
    rq_lengthscale = 1 + 10^(-8), rq_scale = 1
  )
  hp_s <- tibble::tibble(
    rq_variance = 1,
    rq_lengthscale = 1, rq_scale = 1 + 10^(-8)
  )

  deriv_v <- rq_kernel(c(1, 2), c(2, 3), hp, "rq_variance")
  deriv_l <- rq_kernel(c(1, 2), c(2, 3), hp, "rq_lengthscale")
  deriv_s <- rq_kernel(c(1, 2), c(2, 3), hp, "rq_scale")

  emp_deriv_v <- (rq_kernel(c(1, 2), c(2, 3), hp_v)[1] -
                    rq_kernel(c(1, 2), c(2, 3), hp)[1]) / 10^(-8)
  emp_deriv_l <- (rq_kernel(c(1, 2), c(2, 3), hp_l)[1] -
                    rq_kernel(c(1, 2), c(2, 3), hp)[1]) / 10^(-8)
  emp_deriv_s <- (rq_kernel(c(1, 2), c(2, 3), hp_s)[1] -
                    rq_kernel(c(1, 2), c(2, 3), hp)[1]) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
  round(deriv_s, 4) %>% expect_equal(round(emp_deriv_s, 4))
})

test_that("gradients for the Linear kernel are valid", {
  hp <- tibble::tibble(lin_slope = 1, lin_intercept = 1, lin_offset = 1)
  hp_s <- tibble::tibble(
    lin_slope = 1 + 10^(-8),
    lin_intercept = 1, lin_offset = 1
  )
  hp_o <- tibble::tibble(
    lin_slope = 1,
    lin_intercept = 1, lin_offset = 1 + 10^(-8)
  )

  deriv_s <- lin_kernel(c(1, 2), c(2, 3), hp, "lin_slope") %>% as.vector()
  deriv_o <- lin_kernel(c(1, 2), c(2, 3), hp, "lin_offset")

  emp_deriv_s <- (lin_kernel(c(1, 2), c(2, 3), hp_s)[1] -
                    lin_kernel(c(1, 2), c(2, 3), hp)[1]) / 10^(-8)
  emp_deriv_o <- (lin_kernel(c(1, 2), c(2, 3), hp_o)[1] -
                    lin_kernel(c(1, 2), c(2, 3), hp)[1]) / 10^(-8)

  round(deriv_s, 4) %>% expect_equal(round(emp_deriv_s, 4))
  round(deriv_o, 4) %>% expect_equal(round(emp_deriv_o, 4))
})

