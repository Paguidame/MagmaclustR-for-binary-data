
initialisations<- function(prior_mu_k,prior_y_star){
  ## Initialize mu_k according to the value provided by the user
  mu_k <- list()
  if (prior_mu_k %>% is.null()) {
    ## Create a list named by cluster with evaluation of the mean process
    ## at all Input
    for (k in 1:nb_cluster) {
      mu_k[[ID_k[k]]] <- tibble::tibble("Output"=rep(0, length(all_input)),
                                        all_inputs)
    }
    cat(
      "The 'prior_mu_k' argument has not been specified. The hyper_prior mu_k",
      "function is thus set to be 0 everywhere.\n \n"
    )
  } else if (prior_mu_k[[1]] %>% is.function()) {
    ## Create a list named by cluster with evaluation of the mean at all Input
    for (k in 1:nb_cluster) {
      all_inputs %>% dplyr::select(-.data$Reference)
      mu_k[[ID_k[k]]] <- prior_mu_k[[k]](all_inputs)
    }
  } else if (prior_mu_k %>% is.vector()) {
    if (length(prior_mu_k) == nb_cluster) {
      ## Create a list named by cluster with evaluation of the mean at all Input

      for (k in 1:nb_cluster) {
        mu_k[[ID_k[k]]] <- tibble::tibble("Output"=rep(prior_mu_k[[k]],
                                                       length(all_input)),
                                          all_inputs)
      }
    }else if (length(prior_mu_k) == 1) {
      ## Create a list named by cluster with evaluation of the mean at all Input
      for (k in 1:nb_cluster) {
        mu_k[[ID_k[k]]] <- tibble::tibble("Output"=rep(prior_mu_k,
                                                       length(all_input)),
                                          all_inputs)
      }
      cat(
        "The provided 'prior_mu_k' argument is of length 1. Thus, the same",
        "hyper-prior constant mean function has been set for each",
        "cluster.\n \n "
      )
    }else {
      stop(
        "The 'prior_mu_k' argument is of length ", length(prior_mu_k),
        ", whereas there are ", length(hp_k$ID), " clusters."
      )
    }
  }else {
    stop(
      "Incorrect format for the 'prior_mu_k' argument. Please read ",
      "?train_magmaclust() for details."
    )
  }

  # Initialize variables y* for all individuals

  if (is.null(prior_y_star)) {
    all_int <- data %>%
      dplyr::select(-c(.data$ID, .data$Output))
    y_star_prior<- rep(10,length(data$ID))
    y_star<- tibble::tibble("ID"=data$ID,
                            all_int,
                            "Output"=y_star_prior)


  }else if(is.data.frame(prior_y_star)){
    if(!all(names(data) %in% names(prior_y_star))){
      stop("Wrong format for prior_y_star. Make sure that columns are the same
              both in data and prior_y_star.
              Please read ?y_star() for further details.")
    }else {
      y_star<- prior_y_star
    }
  }else{
    stop("The 'prior_y_star' argument must be a data frame. Please read ",
         "?y_star() for further details.")
  }
}

