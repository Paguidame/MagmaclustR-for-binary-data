
#' Training MagmaClust with a Variational EM algorithm
#'
#' The hyper-parameters and the hyper-posterior distributions involved in
#' MagmaClust can be learned thanks to a VEM algorithm implemented in
#' \code{train_magmaclust}. By providing a dataset, the model hypotheses
#' (hyper-prior mean parameters, covariance kernels and number of clusters) and
#' initialisation values for the hyper-parameters, the function computes
#' maximum likelihood estimates of the HPs as well as the mean and covariance
#' parameters of the Gaussian hyper-posterior distributions of the mean
#' processes.
#'
#' @param data A tibble or data frame. Columns required: \code{ID}, \code{Input}
#'    , \code{Output}. Additional columns for covariates can be specified.
#'    The \code{ID} column contains the unique names/codes used to identify each
#'    individual/task (or batch of data).
#'    The \code{Input} column should define the variable that is used as
#'    reference for the observations (e.g. time for longitudinal data). The
#'    \code{Output} column specifies the observed values (the response
#'    variable). The data frame can also provide as many covariates as desired,
#'    with no constraints on the column names. These covariates are additional
#'    inputs (explanatory variables) of the models that are also observed at
#'    each reference \code{Input}.
#' @param nb_cluster A number, indicating the number of clusters of
#'    individuals/tasks that are assumed to exist among the dataset.
#' @param prior_mean_k The set of hyper-prior mean parameters (m_k) for the K
#'    mean GPs, one value for each cluster.
#'     This argument can be specified under various formats, such as:
#'    - NULL (default). All hyper-prior means would be set to 0 everywhere.
#'    - A numerical vector of the same length as the number of clusters.
#'    Each number is associated with one cluster, and considered
#'    to be the hyper-prior mean parameter of the cluster (i.e. a constant
#'    function at all \code{Input}).
#'    - A list of functions. Each function is associated with one cluster. These
#'    functions are all evaluated at all \code{Input} values, to provide
#'    specific hyper-prior mean vectors for each cluster.
#' @param ini_hp_k A tibble or data frame of hyper-parameters
#'    associated with \code{kern_k}, the mean process' kernel.
#'    Required column : \code{ID}. The \code{ID} column contains the unique
#'    names/codes used to identify each cluster. The other columns
#'    should be named according to the hyper-parameters that are used in
#'    \code{kern_k}.
#' @param ini_hp_i A tibble or data frame of hyper-parameters
#'    associated with \code{kern_i}, the individual processes' kernel.
#'    Required column : \code{ID}. The \code{ID} column contains the unique
#'    names/codes used to identify each individual/task. The other columns
#'    should be named according to the hyper-parameters that are used in
#'    \code{kern_i}.
#' @param kern_k A kernel function, associated with the mean GPs.
#'    Several popular kernels
#'    (see \href{https://www.cs.toronto.edu/~duvenaud/cookbook/}{The Kernel
#'    Cookbook}) are already implemented and can be selected within the
#'    following list:
#'    - "SE": (default value) the Squared Exponential Kernel (also called
#'        Radial Basis Function or Gaussian kernel),
#'    - "LIN": the Linear kernel,
#'    - "PERIO": the Periodic kernel,
#'    - "RQ": the Rational Quadratic kernel.
#'    Compound kernels can be created as sums or products of the above kernels.
#'    For combining kernels, simply provide a formula as a character string
#'    where elements are separated by whitespaces (e.g. "SE + PERIO"). As the
#'    elements are treated sequentially from the left to the right, the product
#'    operator '*' shall always be used before the '+' operators (e.g.
#'    'SE * LIN + RQ' is valid whereas 'RQ + SE * LIN' is  not).
#' @param kern_i A kernel function, associated with the individual GPs. (See
#'    details above in \code{kern_k}).
#' @param ini_mixture Initial values of the probability to belong to each
#'    cluster for each individual (\code{\link{ini_mixture}} can be used for
#'    a k-means initialisation. Used by default if NULL).
#' @param common_hp_k A boolean indicating whether hyper-parameters are common
#'    among the mean GPs.
#' @param common_hp_i A boolean indicating whether hyper-parameters are common
#'    among the individual GPs.
#' @param grid_inputs A vector, indicating the grid of additional reference
#'    inputs on which the mean processes' hyper-posteriors should be evaluated.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#' @param n_iter_max A number, indicating the maximum number of iterations of
#'    the VEM algorithm to proceed while not reaching convergence.
#' @param cv_threshold A number, indicating the threshold of the likelihood gain
#'    under which the VEM algorithm will stop. The convergence condition is
#'    defined as the difference of elbo between two consecutive steps,
#'    divided by the absolute value of the last one
#'    ( \eqn{(ELBO_n - ELBO_{n-1}) / |ELBO_n| } ).
#' @param categorial A boolean indicating whether we are in categorial
#'        case for MagmaclustR
#' @param fast_approx A boolean, indicating whether the VEM algorithm should
#'    stop after only one iteration of the VE-step. This advanced feature is
#'    mainly used to provide a faster approximation of the model selection
#'    procedure, by preventing any optimisation over the hyper-parameters.
#'
#' @details The user can specify custom kernel functions for the argument
#'    \code{kern_k} and \code{kern_i}. The hyper-parameters used in the kernel
#'    should have explicit names, and be contained within the \code{hp}
#'    argument. \code{hp} should typically be defined as a named vector or a
#'    data frame. Although it is not mandatory for the \code{train_magmaclust}
#'    function to run, gradients be can provided within kernel function
#'    definition. See for example \code{\link{se_kernel}} to create a custom
#'    kernel function displaying an adequate format to be used in
#'    MagmaClust.
#'
#' @return A list, containing the results of the VEM algorithm used in the
#'    training step of MagmaClust. The elements of the list are:
#'    - hp_k: A tibble containing the trained hyper-parameters for the mean
#'    process' kernel and the mixture proportions for each cluster.
#'    - hp_i: A tibble containing the trained hyper-parameters for the
#'    individual processes' kernels.
#'    - hyperpost: A sub-list containing the parameters of the mean processes'
#'    hyper-posterior distribution, namely:
#'      \itemize{
#'        \item mean: A list of tibbles containing, for each cluster, the
#'              hyper-posterior mean parameters evaluated at each \code{Input}.
#'        \item cov: A list of matrices containing, for each cluster, the
#'              hyper-posterior covariance parameter of the mean process.
#'        \item mixture: A tibble, indicating the mixture probabilities in each
#'              cluster for each individual.
#'      }
#'    - ini_args: A list containing the initial function arguments and values
#'    for the hyper-prior means, the hyper-parameters. In particular, if
#'    those arguments were set to NULL, \code{ini_args} allows us to retrieve
#'    the (randomly chosen) initialisations used during training.
#'    - seq_elbo: A vector, containing the sequence of ELBO values associated
#'    with each iteration.
#'    - converged: A logical value indicated whether the algorithm converged.
#'    - training_time: Total running time of the complete training.
#'
#' @export
#'
#' @examples
#' TRUE
train_magmaclust_cat <- function(data,
                                 nb_cluster = NULL,
                                 prior_mean_k = NULL,
                                 prior_mu_k=NULL,
                                 ini_hp_k = NULL,
                                 ini_hp_i = NULL,
                                 kern_k = "SE",
                                 kern_i = "SE",
                                 ini_mixture = NULL,
                                 common_hp_k = TRUE,
                                 common_hp_i = TRUE,
                                 grid_inputs = NULL,
                                 pen_diag = 1e-10,
                                 n_iter_max = 25,
                                 cv_threshold = 1e-3,
                                 fast_approx = FALSE) {

  ## Check for the correct format of the training data
  if (data %>% is.data.frame()) {
    if (!all(c("ID", "Output", "Input") %in% names(data))) {
      stop(
        "The 'data' argument should be a tibble or a data frame containing ",
        "at least the mandatory column names: 'ID', 'Output' and 'Input'"
      )
    }
  } else {
    stop(
      "The 'data' argument should be a tibble or a data frame containing ",
      "at least the mandatory column names: 'ID', 'Output' and 'Input'"
    )
  }

  ##Convert all non ID columns to double (implicitly throw error if not numeric)
  data = data %>% dplyr::mutate(dplyr::across(- .data$ID, as.double))

  ## Check the number of cluster
  if (nb_cluster %>% is.null()) {
    nb_cluster <- 3
    ID_k <- c("K1", "K2", "K3")
    cat(
      "The number of cluster argument has not been specified. There will",
      "be 3 cluster by default. \n \n"
    )
  }

  ## Retrieve or create the names of the clusters
  if (!is.null(ini_hp_k)) {
    ID_k <- ini_hp_k$ID %>% unique()
    if (length(ID_k) != nb_cluster) {
      stop(
        "The argument 'ini_hp_k' provides hyper-parameters for a number of ",
        "clusters that is different from the 'nb_cluster' argument. "
      )
    }
  } else {
    ID_k <- paste0("K", 1:nb_cluster)
  }

  ## Certify that IDs are of type 'character'
  data$ID <- data$ID %>% as.character()
  ## Extract the list of different IDs
  list_ID <- data$ID %>% unique()

  ## Get input column names
  names_col <- data %>%
    dplyr::select(-c(.data$ID,.data$Output)) %>%
    names()

  ## Keep 6 significant digits for entries to avoid numerical errors and
  ## Add a Reference column for identification and sort according to it
  data <- data %>% purrr::modify_at(tidyselect::all_of(names_col),signif) %>%
    tidyr::unite("Reference",
                 tidyselect::all_of(names_col),
                 sep=":",
                 remove = FALSE) %>%
    tidyr::drop_na() %>%
    dplyr::group_by(.data$ID) %>%
    dplyr::arrange(.data$Reference, .by_group = TRUE) %>%
    dplyr::ungroup()

  ## Check that individuals do not have duplicate inputs
  if(!(setequal(data %>% dplyr::select(-.data$Output),
                data %>% dplyr::select(-.data$Output) %>% unique() ))
  ){
    stop("At least one individual have several Outputs on the same grid point.",
         " Please read ?train_magma() for further details."
    )
  }

  ## Extract the union of all reference inputs provided in the training data
  all_input <- data %>%
    dplyr::pull(.data$Reference) %>%
    unique() %>%
    sort()

  all_inputs <- data %>%
    dplyr::select(-c(.data$ID, .data$Output)) %>%
    unique() %>%
    dplyr::arrange(.data$Reference)

  ## Initialise the individual process' HPs according to user's values
  if (kern_i %>% is.function()) {
    if (ini_hp_i %>% is.null()) {
      stop(
        "When using a custom kernel function the 'ini_hp_i' argument is ",
        "mandatory, in order to provide the name of the hyper-parameters. ",
        "You can use the function 'hp()' to easily generate a tibble of random",
        " hyper-parameters with the desired format for initialisation."
      )
    }
  } else {
    if (ini_hp_i %>% is.null()) {
      hp_i <- hp(kern_i,
                 list_ID = list_ID,
                 common_hp = common_hp_i,
                 noise = TRUE
      )
      cat(
        "The 'ini_hp_i' argument has not been specified. Random values of",
        "hyper-parameters for the individual processes are used as",
        "initialisation.\n \n"
      )
    } else if (!("ID" %in% names(ini_hp_i))) {
      ## Create a full tibble of common HPs if the column ID is not specified
      hp_i <- tibble::tibble(
        ID = list_ID,
        dplyr::bind_rows(ini_hp_i)
      )
      cat(
        "No 'ID' column in the 'ini_hp_i' argument. The same hyper-parameter",
        "values have been duplicated for every 'ID' present in the 'data'.\n \n"
      )
    } else if (!(all(as.character(ini_hp_i$ID) %in% as.character(list_ID)) &
                 all(as.character(list_ID) %in% as.character(ini_hp_i$ID)))) {
      stop(
        "The 'ID' column in 'ini_hp_i' is different from the 'ID' of the ",
        "'data'."
      )
    } else {
      hp_i <- ini_hp_i
    }
  }

  ## Add a 'noise' hyper-parameter if absent
  if (!("noise" %in% names(hp_i))) {
    if (common_hp_i) {
      hp_i <- hp_i %>% dplyr::mutate(hp(NULL, noise = T))
    }else {
      hp_i <- hp_i %>%
        dplyr::left_join(hp(NULL,
                            list_ID = hp_i$ID,
                            noise = T),
                         by = "ID"
        )
    }
  }

  ## Initialise the cluster process' hp according to user's values
  if (kern_k %>% is.function()) {
    if (ini_hp_k %>% is.null()) {
      stop(
        "When using a custom kernel function the 'ini_hp_k' argument is ",
        "mandatory, in order to provide the name of the hyper-parameters. ",
        "You can use the function 'hp()' to easily generate a tibble of random",
        " hyper-parameters with the desired format for initialisation."
      )
    }
  } else {
    if (ini_hp_k %>% is.null()) {
      hp_k <- hp(kern_k,
                 list_ID = ID_k,
                 common_hp = common_hp_k,
                 noise = F
      )
      cat(
        "The 'ini_hp_k' argument has not been specified. Random values of",
        "hyper-parameters for the mean processes are used as",
        "initialisation.\n \n"
      )
    } else if (!("ID" %in% names(ini_hp_k))) {
      ## Create a full tibble of common HPs if the column ID is not specified
      hp_k <- tibble::tibble(
        'ID' = ID_k,
        dplyr::bind_rows(ini_hp_k)
      )
      cat(
        "No 'ID' column in the 'ini_hp_k' argument. The same hyper-parameter",
        "values have been duplicated for every cluster's 'ID'.\n \n"
      )
    } else {
      hp_k <- ini_hp_k
    }
  }

  ## Initialise m_k according to the value provided by the user
  m_k <- list()
  if (prior_mean_k %>% is.null()) {
    ## Create a list named by cluster with evaluation of the mean at all Input
    for (k in 1:nb_cluster) {
      m_k[[ID_k[k]]] <- rep(0, length(all_input))
    }
    cat(
      "The 'prior_mean' argument has not been specified. The hyper_prior mean",
      "function is thus set to be 0 everywhere.\n \n"
    )
  } else if (prior_mean_k[[1]] %>% is.function()) {
    ## Create a list named by cluster with evaluation of the mean at all Input
    for (k in 1:nb_cluster) {
      all_inputs %>% dplyr::select(-.data$Reference)
      m_k[[ID_k[k]]] <- prior_mean_k[[k]](all_inputs)
    }
  } else if (prior_mean_k %>% is.vector()) {
    if (length(prior_mean_k) == nb_cluster) {
      ## Create a list named by cluster with evaluation of the mean at all Input
      for (k in 1:nb_cluster) {
        m_k[[ID_k[k]]] <- rep(prior_mean_k[[k]], length(all_input))
      }
    } else if (length(prior_mean_k) == 1) {
      ## Create a list named by cluster with evaluation of the mean at all Input
      for (k in 1:nb_cluster) {
        m_k[[ID_k[k]]] <- rep(prior_mean_k, length(all_input))
      }
      cat(
        "The provided 'prior_mean' argument is of length 1. Thus, the same",
        "hyper-prior constant mean function has been set for each",
        "cluster.\n \n "
      )
    }else {
      stop(
        "The 'prior_mean_k' argument is of length ", length(prior_mean_k),
        ", whereas there are ", length(hp_k$ID), " clusters."
      )
    }
  }else {
    stop(
      "Incorrect format for the 'prior_mean_k' argument. Please read ",
      "?train_magmaclust() for details."
    )
  }

  ## Initialize mu_k according to the value provided by the user
  mu_k <- list()
  if (prior_mu_k %>% is.null()) {
    ## Create a list named by cluster with evaluation of the mean process
    ## at all Input
    for (k in 1:nb_cluster) {
      ## mu_k[[ID_k[k]]] <- tibble::tibble("Output"= rmvnorm(1, mean =   rep(0, length(all_input)),  diag(length(all_input)))%>% as.vector(),
      ##                               all_inputs)
      #mu_k[[ID_k[k]]] <- tibble::tibble("Output"= base::sample(c(-1,0, 1), length(all_input), replace = TRUE), all_inputs)
      mu_k[[ID_k[k]]] <- tibble::tibble("Output"= rep(0, length(all_input)), all_inputs)
    }
    cat(
      "The 'prior_mu_k' argument has not been specified. The hyper_prior mu_k",
      "function is thus set to be 0 everywhere.\n \n"
    )
  }else if (prior_mu_k[[1]] %>% is.function()) {
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

  ## Track the total training time
  t1 <- Sys.time()

  ## Initialize the monitoring information
  cv <- FALSE
  elbo_monitoring <- -Inf
  seq_elbo <- c()

  if (is.null(ini_mixture)) {
    mixture <- ini_mixture(data,
                           k = nb_cluster,
                           name_clust = ID_k,
                           50)
  }else if(is.data.frame(ini_mixture)){
    if(!all(c("ID", ID_k) %in% names(ini_mixture))){
      stop("Wrong format for ini_mixture. Make sure that the number of ",
           "clusters are the same both in 'train_magmaclust()' and ",
           "ini_mixture. Please read ?ini_mixture() for further details.")
    }else {
      mixture <- ini_mixture
    }
  }else{
    stop("The 'ini_mixture' argument must be a data frame. Please read ",
         "?ini_mixture() for further details.")
  }


  # Initialize the affectations variables

  old_z <- simu_affectations(mixture)
  #while (sum(old_z[,"K1"]) < 2 |sum(old_z[,"K2"]) < 2 | sum(old_z[,"K3"]) < 2){
  # old_z <- simu_affectations(mixture)
  #}

  hp_k[["prop_mixture"]] <- old_z %>%
    dplyr::select(-.data$ID) %>%
    colMeans() %>%
    as.vector()

  ## Keep an history of the (possibly random) initial values of hyper-parameters
  hp_i_ini <- hp_i
  hp_k_ini <- hp_k
  mixture_ini <- mixture
  ## Iterate VE-step and VM-step until convergence
  for (i in 1:n_iter_max)
  {
    ## Track the running time for each iteration of the EM algorithm
    t_i_1 <- Sys.time()
    # run the SEM algorithm in case of categorial for MagmaclustR

    ## SE-Step of MagmaClust
    post <- se_step2(
      data,
      mu_k,
      m_k,
      kern_k,
      kern_i,
      hp_k,
      hp_i,
      old_z,
      iter = i,
      pen_diag
    )
    ## Break after VE-step if we can to compute the fast approximation
    if (fast_approx) {
      ## Track the ELBO values
      seq_elbo <- elbo_monitoring_VEM(
        hp_k,
        hp_i,
        data,
        kern_i,
        kern_k,
        hyperpost=NULL,
        m_k = m_k,
        latents = post[c( "Z","y_star","mu")],
        pen_diag,
        categorial = TRUE
      )

      cv <- FALSE
      break
    }

    ## SM-Step of MagmaClsut
    new_hp <- sm_step2(post$y_star,
                       hp_k,
                       hp_i,
                       list_latents = list("mu"= post$mu,"Z"= post$Z),
                       kern_k,
                       kern_i,
                       m_k,
                       common_hp_k,
                       common_hp_i,
                       pen_diag
    ) # %>% suppressMessages()

    new_hp_k <- new_hp$hp_k
    new_hp_i <- new_hp$hp_i

    ## In case something went wrong during the optimization
    if (any(is.na(new_hp_k)) | any(is.na(new_hp_i))) {
      warning(paste0("The M-step encountered an error at iteration : ", i))
      warning(
        "Training has stopped and hyper-parameters values from the ",
        "last valid iteration are returned."
      )
      break
    }

    ## Monitoring of the elbo
    new_elbo_monitoring <- elbo_monitoring_VEM(
      new_hp_k,
      new_hp_i,
      data,
      kern_i,
      kern_k,
      hyperpost=NULL,
      m_k = m_k,
      pen_diag,
      latents = post[c( "Z","y_star","mu")],
      categorial = TRUE
    )

    diff_moni <- new_elbo_monitoring - elbo_monitoring
    if (diff_moni %>% is.nan()) {
      diff_moni <- -Inf
    }

    #if (diff_moni < 0) {
    #warning("Likelihood descreased")
    #}

    ## Update HPs values and the elbo monitoring
    hp_k <- new_hp_k
    hp_i <- new_hp_i
    mu_k <- post$mu
    #mixture <- post$mixture
    old_z <- post$Z
    elbo_monitoring <- new_elbo_monitoring

    ## Track the ELBO values
    seq_elbo <- c(seq_elbo, elbo_monitoring)

    ## Compute the convergence ratio
    eps <- diff_moni / abs(elbo_monitoring)
    if (eps %>% is.nan()) {
      eps <- 1
    }

    ## Provide monitoring information
    t_i_2 <- Sys.time()
    paste0(
      "SEM algorithm, step ", i, ": ",
      difftime(t_i_2, t_i_1, units = "secs") %>% round(2),
      " seconds \n \n"
    ) %>%
      cat()

    paste0(
      "Value of the elbo: ",
      elbo_monitoring %>% round(5),
      " --- Convergence ratio = ",
      eps %>% round(5),
      "\n \n"
    ) %>%
      cat()

    ## Check the convergence condition

    if (abs(eps) < cv_threshold) {
      cat(
        "The SEM algorithm successfully converged, training is completed.",
        "\n \n"
      )
      cv <- TRUE
      break
    }
    ## Check for a prematurate ending of the EM algorithm
    if (!cv & (i == n_iter_max)) {
      warning(
        "The SEM algorithm has reached the maximum number of iterations ",
        "before convergence, training might be sub-optimal \n \n"
      )
    }
  }

  ## Evaluate the hyper-posterior on the grid of inputs if provided
  if (!is.null(grid_inputs)) {
    cat(
      "Start evaluating hyper-posterior distributions of the mean processes",
      "on the provided grid of inputs... \n \n"
    )

    post <- hyperposterior_clust(
      data = data,
      post$mixture,
      hp_k = hp_k,
      hp_i = hp_i,
      kern_k = kern_k,
      kern_i = kern_i,
      prior_mean_k = prior_mean_k,
      grid_inputs = grid_inputs,
      pen_diag = pen_diag
    )
    cat("Done!\n \n")
  } else {
    ## Create a variable for directly plotting the mean process' hyper-posterior
    floop_pred <- function(k) {
      tibble::tibble(post$mean[[k]],
                     "Var" =  post$cov[[k]] %>%
                       diag() %>%
                       as.vector()
      ) %>%
        dplyr::rename("Mean" = .data$Output) %>%
        dplyr::select(-.data$Reference) %>%
        return()
    }
    post$pred <- sapply(ID_k, floop_pred, simplify = FALSE, USE.NAMES = TRUE)
  }


  ## Create an history list of the initial arguments of the function
  fct_args <- list(
    "data" = data %>% dplyr::select(-.data$Reference),
    "nb_cluster" = nb_cluster,
    "prior_mean_k" = prior_mean_k,
    "ini_hp_k" = hp_k_ini,
    "ini_hp_i" = hp_i_ini,
    "kern_k" = kern_k,
    "kern_i" = kern_i,
    "ini_mixture" = mixture_ini,
    "common_hp_k" = common_hp_k,
    "common_hp_i" = common_hp_i,
    "n_iter_max" = n_iter_max,
    "pen_diag" = pen_diag,
    "cv_threshold" = cv_threshold
  )
  t2 <- Sys.time()
  ## Create and return the list of elements from the trained model
  list(
    "hp_k" = hp_k,
    "hp_i" = hp_i,
    "hyperpost" = post,
    "ini_args" = fct_args,
    "seq_elbo" = seq_elbo,
    "converged" = cv,
    "training_time" = difftime(t2, t1, units = "secs")
  ) %>%
    return()
}

