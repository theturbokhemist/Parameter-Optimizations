##############################generate_ref_object########################################################################################################################################################
generate_ref_object <- function(topology, reference_state, numModels = 10000, integrateStepSize = 0.02, global_parameter) {
  
  reference_rset <- sracipeSimulate(topology, integrate = FALSE,
                                    numModels = numModels, genParams = TRUE,
                                    integrateStepSize = integrateStepSize)
  
  state <- reference_state
  
  #hill gene indexes
  hill_genes <- grep(pattern = "N", colnames(sracipeParams(reference_rset)))
  
  if (!missing(global_parameter)) {
    
    parametric_variation_reference_state <- reference_state
    
    for (w in 1:length(parametric_variation_reference_state)) {
      
      new_bounds <- parametric_variation_bounds(parameter_bounds = reference_state[[w]], global_parameter = global_parameter)
      
      parametric_variation_reference_state[[w]][[1]] <- new_bounds$bound_min
      parametric_variation_reference_state[[w]][[2]] <- new_bounds$bound_max
      
    }
    
    for (u in 1:length(hill_genes)) {
      
      parametric_variation_reference_state[[hill_genes[u]]][[1]] <- round(parametric_variation_reference_state[[hill_genes[u]]][[1]])
      parametric_variation_reference_state[[hill_genes[u]]][[2]] <- round(parametric_variation_reference_state[[hill_genes[u]]][[2]]) 
      
    }
    
    state <- parametric_variation_reference_state
    
  }
  
  
  #assign vaues from uniform samples distributions to each model of reference rset.  
  for (i in 1:length(state)) {
    
    sracipeParams(reference_rset)[,i] <- runif(numModels, min= as.numeric(state[[i]][1]), max= as.numeric(state[[i]][2]))  
    
  }
  
  for (k in 1:length(hill_genes)) {
    #you need to subtract 0.5 from lower bound and add 0.5 to upper bound in order to get a uniform distribution
    sracipeParams(reference_rset)[,hill_genes[k]] <- round(runif(numModels, min= as.numeric(state[[hill_genes[k]]][1]) - 0.5, max= as.numeric(state[[hill_genes[k]]][2]) + 0.5))  
    
  }
  
  ref_expression_rset <- sracipeSimulate(reference_rset, integrate = TRUE, genParams = FALSE)
  
  expr_reference_log <- log2(t(assay(ref_expression_rset)))
  
  means_reference <- colMeans(expr_reference_log)
  
  sds_reference <- apply(expr_reference_log, 2, sd)
  
  #z normalize
  expr_reference_log_z <- sweep(expr_reference_log, 
                                2, means_reference, FUN = "-")
  
  expr_reference_log_z <- sweep(expr_reference_log_z, 
                                2, sds_reference, FUN = "/")
  
  #find eigenvector,loading scores
  prcomp_reference <- prcomp(expr_reference_log_z, center = FALSE, scale. = FALSE)
  
  #list of all the ref information
  reference_rset_object <- list(reference_rset = ref_expression_rset, reference_state = reference_state, expression_ref = expr_reference_log_z, prcomp_ref = prcomp_reference, means_ref = means_reference, sds_ref = sds_reference, topology = topology)
  
  if (!missing(global_parameter)) {
    
    reference_rset_object <- list(reference_rset = ref_expression_rset, reference_state = reference_state, expression_ref = expr_reference_log_z, prcomp_ref = prcomp_reference, means_ref = means_reference, sds_ref = sds_reference, topology = topology,
                                  parametric_variation_reference_state = parametric_variation_reference_state, global_parameter = global_parameter)
    
  }
  
  reference_rset_object
  
}

#####################calculate_score##################################################################################################################################################################

calc_score <- function(state, reference_rset_object, scoring_method = "Default", numModels = 10000, integrateStepSize = 0.02, save_rset = FALSE, filter) {
  
  parameter_names <- colnames(sracipeParams(reference_rset_object$reference_rset))
  
  hill_params <- grep(pattern = "N", parameter_names)
  
  simulated_rset <- sracipeSimulate(reference_rset_object$topology, integrate = FALSE,
                                    numModels = numModels, genParams = TRUE,
                                    integrateStepSize = integrateStepSize)
  
  for (t in 1:length(state)) {
    
    sracipeParams(simulated_rset)[,t] <- runif(numModels, min = state[[t]][[1]], max = state[[t]][[2]])  
    
  }
  
  #sample and round all the hill coefficients
  for (u in 1:length(hill_params)) {
    #you need to subtract 0.5 from lower bound and add 0.5 to upper bound in order to get a uniform distribution
    sracipeParams(simulated_rset)[,hill_params[u]] <- round(runif(numModels, min= round(state[[hill_params[u]]][[1]]) - 0.5, max= round(state[[hill_params[u]]][[2]]) + 0.5))  
    
  }
  
  simulated_rset <- sracipeSimulate(simulated_rset, integrate = TRUE, genParams = FALSE)
  
  simulated_rset_expr_log <- log2(t(assay(simulated_rset)))
  
  simulated_rset_expr_log_z <- sweep(simulated_rset_expr_log,
                                     2, reference_rset_object$means_ref, FUN = "-")
  
  simulated_rset_expr_log_z <- sweep(simulated_rset_expr_log_z,
                                     2, reference_rset_object$sds_ref, FUN = "/")
  
  #rotate
  simulated_rset_expr_PCs <- simulated_rset_expr_log_z %*% reference_rset_object$prcomp_ref$rotation
  
  #cost_function vector
  cf_vec <- vector(length = ncol(simulated_rset_expr_log_z))
  
  #scoring
  if (scoring_method1 == "PCs") {
    
    simulated_expression <- simulated_rset_expr_PCs
    
    #ks score for distribution of each gene
    for (w in 1:length(cf_vec)) {
      
      ks_test <- dgof::ks.test(x = reference_rset_object$prcomp_ref$x[,w], y = simulated_rset_expr_PCs[,w])
      cf_vec[w] <- as.numeric(ks_test$statistic)
      
    }
    
  } else {
    
    simulated_expression <- simulated_rset_expr_log_z
    
    for (w in 1:length(cf_vec)) {
      
      ks_test <- dgof::ks.test(x = reference_rset_object$expression_ref[,w], y = simulated_rset_expr_log_z[,w])
      cf_vec[w] <- as.numeric(ks_test$statistic)
      
    }
  }
  
  costs_object <- list(cf_vec = cf_vec, cf_vec_sum = sum(cf_vec), simulated_rset = simulated_rset)
  
  if (save_rset == TRUE) {
    
    costs_object <- list(cf_vec = cf_vec, cf_vec_sum = sum(cf_vec), simulated_rset = simulated_rset,  simulated_expression = simulated_expression,
                         arguments = list(scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, state = state))
    
  }
  costs_object
}



#####################filter_score##################################################################################################################################################################
filter_score <- function(current_state_rset_object, reference_rset_object, scoring_method = "Default", filter) {

  
  
  filtered_expr_log <- log2(t(assay(current_state_rset_object$simulated_rset)))[filter, ]
  
  filtered_expr_log_z <- sweep(filtered_expr_log,
                               2, reference_rset_object$means_ref, FUN = "-")
  
  filtered_expr_log_z <- sweep(filtered_expr_log_z,
                               2, reference_rset_object$sds_ref, FUN = "/")
  
  #rotate
  filtered_expr_PCs <- filtered_expr_log_z %*% reference_rset_object$prcomp_ref$rotation
  
  #cost_function vector
  cf_vec <- vector(length = ncol(filtered_expr_log_z))
  
  #scoring
  if (scoring_method == "PCs") {
    
    filtered_expression <- filtered_expr_PCs
    
    #ks score for distribution of each gene
    for (w in 1:length(cf_vec)) {
      
      ks_test <- dgof::ks.test(x = reference_rset_object$prcomp_ref$x[,w], y = filtered_expr_PCs[,w])
      cf_vec[w] <- as.numeric(ks_test$statistic)
      
    }
    
  } else {
    
    filtered_expression <- filtered_expr_log_z
    
    for (w in 1:length(cf_vec)) {
      
      ks_test <- dgof::ks.test(x = reference_rset_object$expression_ref[,w], y = filtered_expr_log_z[,w])
      cf_vec[w] <- as.numeric(ks_test$statistic)
      
    }
  }
  
  costs_object <- list(cf_vec = cf_vec, cf_vec_sum = sum(cf_vec), filtered_expression = filtered_expression)
  

  costs_object
  
}


#####################simulate_partial_state##################################################################################################################################################################
simulate_partial_state <- function(partial_state, current_state_rset_object, reference_rset_object, partial_rset, scoring_method = "Default", hill_params, param, bound, save_all = FALSE) {
  
  numModels <- nrow(sracipeParams(partial_rset))
  
  #assign vaues from uniform samples distributions to each model of reference rset.  
  for (i in 1:length(partial_state)) {
    
    sracipeParams(partial_rset)[,i] <- runif(numModels, min = partial_state[[i]][[1]], max = partial_state[[i]][[2]])  
    
  }
  
  #replace values for specific parameter
  if (bound == 1) {
    
    sracipeParams(partial_rset)[,param] <- runif(numModels, min = partial_state[[param]][[1]], max = current_state_rset_object$arguments$state[[param]][[1]])
    
  } else {
    
    sracipeParams(partial_rset)[,param] <- runif(numModels, min = current_state_rset_object$arguments$state[[param]][[2]], max = partial_state[[param]][[2]])
    
    
  }
  
  for (k in 1:length(hill_params)) {
    #you need to subtract 0.5 from lower bound and add 0.5 to upper bound in order to get a uniform distribution
    sracipeParams(partial_rset)[,hill_params[k]] <- round(runif(numModels, min = round(partial_state[[hill_params[k]]][[1]]) - 0.5, max = round(partial_state[[hill_params[k]]][[2]]) + 0.5))  
    
  }
  
  partial_state_expression_rset <- sracipeSimulate(partial_rset, integrate = TRUE, genParams = FALSE)
  
  partial_state_expr_log <- rbind(log2(t(assay(current_state_rset_object$simulated_rset))), log2(t(assay(partial_state_expression_rset))))
  
  partial_state_expr_log_z <- sweep(partial_state_expr_log,
                               2, reference_rset_object$means_ref, FUN = "-")
  
  partial_state_expr_log_z <- sweep(partial_state_expr_log_z,
                               2, reference_rset_object$sds_ref, FUN = "/")
  
  #rotate
  partial_state_expr_PCs <- partial_state_expr_log_z %*% reference_rset_object$prcomp_ref$rotation
  
  #cost_function vector
  cf_vec <- vector(length = ncol(partial_state_expr_log_z))
  
  #scoring
  if (scoring_method == "PCs") {
    
    partial_state_expression <- partial_state_expr_PCs
    
    #ks score for distribution of each gene
    for (w in 1:length(cf_vec)) {
      
      ks_test <- dgof::ks.test(x = reference_rset_object$prcomp_ref$x[,w], y = partial_state_expr_PCs[,w])
      cf_vec[w] <- as.numeric(ks_test$statistic)
      
    }
    
  } else {
    
    partial_state_expression <- partial_state_expr_log_z
    
    for (w in 1:length(cf_vec)) {
      
      ks_test <- dgof::ks.test(x = reference_rset_object$expression_ref[,w], y = partial_state_expr_log_z[,w])
      cf_vec[w] <- as.numeric(ks_test$statistic)
      
    }
  }
  
  costs_object <- list(cf_vec = cf_vec, cf_vec_sum = sum(cf_vec), partial_state_expression = partial_state_expression)
  
  if (save_all == TRUE) {
    
    costs_object <- list(cf_vec = cf_vec, cf_vec_sum = sum(cf_vec), partial_state_expression = partial_state_expression, partial_state_expression_rset = partial_state_expression_rset)
    
  }
  
  costs_object

}


#####################calc_partial_derivs##################################################################################################################################################################

calc_partial_derivatives <- function(current_state, reference_rset_object, global_bounds, step_size = 0.01, numModels_filter_state = 20000, integrateStepSize_filter_state = 0.02, scoring_method = scoring_method, save_all = FALSE) {
  
  #number of parameters
  number_of_params <- length(current_state)
  
  #vector of paramater names
  parameter_names <- colnames(sracipeParams(reference_rset_object$reference_rset))
  
  #vector of parameter_names + bound
  paramater_bounds_names <- paste(rep(parameter_names, each = 2),
                            rep(c("lower", "upper"), times = number_of_params), sep = "_")
  
  #index of each parameter type in vector of parameter names
  production_params <- grep(pattern = "G", parameter_names)
  degradation_params <- grep(pattern = "K", parameter_names)
  threshold_params <- grep(pattern = "TH", parameter_names)
  hill_params <- grep(pattern = "N", parameter_names)
  FC_params <- grep(pattern = "FC", parameter_names)
  
  #vector of step magnitudes for each parameter. The step magnitude is simply the step_size *( upper bound -  lowerbound) 
  update_vector <- vector(length = number_of_params)
  update_vector_rep <- vector(length = (number_of_params)*4)
  
  counter <- 0
    
    for (i in 1:number_of_params) {
      
      update_vector[i] <- (current_state[[i]][[2]] - current_state[[i]][[1]])*step_size
      
      update_vector_rep[(i + counter):(i+counter+3)] <- rep(update_vector[i], 4)
      
      counter <- counter + 3
      
    }
  
  #list of the current states bounds for each parameter *4
  current_state_bounds_list <- vector(mode = "list", length = (number_of_params)*4)
  
  counter <- 0
  i <- 1
  
  for (i in 1:number_of_params) {
    
    current_state_bounds_list[[i + counter]] <- current_state[[i]]
    current_state_bounds_list[[i + counter + 1]] <- current_state[[i]]
    current_state_bounds_list[[i + counter + 2]] <- current_state[[i]]
    current_state_bounds_list[[i + counter + 3]] <- current_state[[i]]
    
    counter <- counter + 3
    
  }

  #generate current_state_score_object
  current_state_score_object <- calc_score(state = current_state, reference_rset_object = reference_rset_object, scoring_method = scoring_method, numModels = numModels_filter_state, integrateStepSize = integrateStepSize_filter_state, save_rset = TRUE)
  current_state_params <- sracipeParams(current_state_score_object$simulated_rset)
  
  #generate list of partial states
  partial_states_list <- vector(mode = "list", length = number_of_params*4) #*2 for 2 bounds and *2 for adding and subtracting from each bound
  
  counter <- 0
  
  for (j in 1:length(current_state)) {
    
    #first lower/upper refers to which bound, second minus/plus refers t
    
    partial_state_lower_minus <- current_state
    partial_state_lower_plus <- current_state
    
    partial_state_upper_minus <- current_state
    partial_state_upper_plus <- current_state
    
    partial_state_lower_minus[[j]][1] <- current_state[[j]][[1]] - update_vector[j]
    partial_state_lower_plus[[j]][1] <- current_state[[j]][[1]] + update_vector[j]
    
    partial_state_upper_minus[[j]][2] <- current_state[[j]][[2]] - update_vector[j]
    partial_state_upper_plus[[j]][2] <- current_state[[j]][[2]] + update_vector[j]
    
    
    partial_states_list[[j + counter]] <- partial_state_lower_minus
    partial_states_list[[j + counter + 1]] <- partial_state_lower_plus
    partial_states_list[[j + counter + 2]] <- partial_state_upper_minus
    partial_states_list[[j + counter + 3]] <- partial_state_upper_plus
    
    
    counter <- counter + 3
    
  }
  
  
  
  counter <- 0

  partial_derivs_vec <- vector(length = length(partial_states_list))
  percent_models <- vector(length = length(partial_states_list))
  
  filters_list <- vector(mode = "list", length = length(partial_states_list))
  
  #generate parameter for amount of models that will be added to current state if partial state bounds are outside of current state bounds
  partial_rset <- sracipeSimulate(reference_rset_object$topology, integrate = FALSE,
                                  numModels = round(numModels_filter_state*step_size), genParams = TRUE,
                                  integrateStepSize = integrateStepSize_filter_state)
  
  
  for (k in 1:length(current_state)) {
    
    if ((partial_states_list[[k + counter]][[k]][[1]] >= current_state[[k]][[1]]) & (partial_states_list[[k + counter]][[k]][[2]] <= current_state[[k]][[2]])) {
      
      #filter models of the current state using the param bounds of the partial state
      
      filtered_models <- which(current_state_params[,k] >= partial_states_list[[k + counter]][[k]][[1]] & current_state_params[,k] <= partial_states_list[[k + counter]][[k]][[2]])
      
      #round the partial states if the hill coefficient was changed
      if (k %in% hill_params) {
        
        filtered_models <- which(current_state_params[,k] >= round(partial_states_list[[k + counter]][[k]][[1]]) & current_state_params[,k] <= round(partial_states_list[[k + counter]][[k]][[2]]))
      }
      
      filtered_scores_object <- filter_score(current_state_rset_object = current_state_score_object, reference_rset_object = reference_rset_object, scoring_method = scoring_method, filter = filtered_models)
      
      if ((partial_states_list[[k + counter]][[k]][[2]] < current_state[[k]][[2]])) {
        
        partial_deriv <- -(filtered_scores_object$cf_vec_sum - current_state_score_object$cf_vec_sum)/update_vector_rep[k + counter]
        
      } else {
        
        partial_deriv <- (filtered_scores_object$cf_vec_sum - current_state_score_object$cf_vec_sum)/update_vector_rep[k + counter]
        
      }
      
      partial_derivs_vec[k + counter] <- partial_deriv
      
      percent_models[k + counter] <- 100*length(filtered_models)/numModels_filter_state
      
      filters_list[[k + counter]] <- filtered_models
      
    } else {
      
      if ((partial_states_list[[k + counter]][[k]][[1]] < current_state[[k]][[1]])) {
        
        bound <- 1
        
      } else {
        
        bound <- 2
        
      }
      
      partial_state_costs_object <- simulate_partial_state(partial_state = partial_states_list[[k + counter]], current_state_rset_object = current_state_score_object,
                                                           reference_rset_object = reference_rset_object, partial_rset = partial_rset, scoring_method = scoring_method, hill_params = hill_params, param = k, bound = bound)
      
      if (bound == 1) {
        
        partial_deriv <- -(partial_state_costs_object$cf_vec_sum - current_state_score_object$cf_vec_sum)/update_vector_rep[k + counter]
        
      } else {
        
        partial_deriv <- (partial_state_costs_object$cf_vec_sum - current_state_score_object$cf_vec_sum)/update_vector_rep[k + counter]
        
      }
      
      
      partial_derivs_vec[k + counter] <- partial_deriv
      percent_models[k + counter] <- 100*nrow(partial_state_costs_object$partial_state_expression)/numModels_filter_state
      
    }
    
    ########################################
    if ((partial_states_list[[k + counter + 1]][[k]][[1]] >= current_state[[k]][[1]]) & (partial_states_list[[k + counter + 1]][[k]][[2]] <= current_state[[k]][[2]])) {
      
      #filter models of the current state using the param bounds of the partial state

      filtered_models <- which(current_state_params[,k] >= partial_states_list[[k + counter + 1]][[k]][[1]] & current_state_params[,k] <= partial_states_list[[k + counter + 1]][[k]][[2]])
      
      #round the partial states if the hill coefficient was changed
      if (k %in% hill_params) {
        
        filtered_models <- which(current_state_params[,k] >= round(partial_states_list[[k + counter + 1]][[k]][[1]]) & current_state_params[,k] <= round(partial_states_list[[k + counter + 1]][[k]][[2]]))
      }
      
      filtered_scores_object <- filter_score(current_state_rset_object = current_state_score_object, reference_rset_object = reference_rset_object, scoring_method = scoring_method, filter = filtered_models)
      
      if ((partial_states_list[[k + counter + 1]][[k]][[2]] < current_state[[k]][[2]])) {
        
        partial_deriv <- -(filtered_scores_object$cf_vec_sum - current_state_score_object$cf_vec_sum)/update_vector_rep[k + counter + 1]
        
      } else {
        
        partial_deriv <- (filtered_scores_object$cf_vec_sum - current_state_score_object$cf_vec_sum)/update_vector_rep[k + counter + 1]
        
      }
      
      partial_derivs_vec[[k + counter + 1]] <- partial_deriv
      percent_models[k + counter + 1] <- 100*length(filtered_models)/numModels_filter_state
      
      filters_list[[k + counter + 1]] <- filtered_models
      
    } else {
      
      if ((partial_states_list[[k + counter + 1]][[k]][[1]] < current_state[[k]][[1]])) {
        
        bound <- 1
        
      } else {
        
        bound <- 2
        
      }
      
      partial_state_costs_object <- simulate_partial_state(partial_state = partial_states_list[[k + counter + 1]], current_state_rset_object = current_state_score_object,
                                                           reference_rset_object = reference_rset_object, partial_rset = partial_rset, scoring_method = scoring_method, hill_params = hill_params, param = k, bound = bound)
      
      if (bound == 1) {
        
        partial_deriv <- -(partial_state_costs_object$cf_vec_sum - current_state_score_object$cf_vec_sum)/update_vector_rep[k + counter + 1]
        
      } else {
        
        partial_deriv <- (partial_state_costs_object$cf_vec_sum - current_state_score_object$cf_vec_sum)/update_vector_rep[k + counter + 1]
        
      }
      
      partial_derivs_vec[k + counter + 1] <- partial_deriv
      percent_models[k + counter + 1] <- 100*nrow(partial_state_costs_object$partial_state_expression)/numModels_filter_state
      
    }
    ########################################
    if ((partial_states_list[[k + counter + 2]][[k]][[1]] >= current_state[[k]][[1]]) & (partial_states_list[[k + counter + 2]][[k]][[2]] <= current_state[[k]][[2]])) {
      
      #filter models of the current state using the param bounds of the partial state
      
      filtered_models <- which(current_state_params[,k] >= partial_states_list[[k + counter + 2]][[k]][[1]] & current_state_params[,k] <= partial_states_list[[k + counter + 2]][[k]][[2]])
      
      #round the partial states if the hill coefficient was changed
      if (k %in% hill_params) {
        
        filtered_models <- which(current_state_params[,k] >= round(partial_states_list[[k + counter + 2]][[k]][[1]]) & current_state_params[,k] <= round(partial_states_list[[k + counter + 2]][[k]][[2]]))
      }
      
      filtered_scores_object <- filter_score(current_state_rset_object = current_state_score_object, reference_rset_object = reference_rset_object, scoring_method = scoring_method, filter = filtered_models)
      
      if ((partial_states_list[[k + counter + 2]][[k]][[2]] < current_state[[k]][[2]])) {
        
        partial_deriv <- -(filtered_scores_object$cf_vec_sum - current_state_score_object$cf_vec_sum)/update_vector_rep[k + counter + 2]
        
      } else {
        
        partial_deriv <- (filtered_scores_object$cf_vec_sum - current_state_score_object$cf_vec_sum)/update_vector_rep[k + counter + 2]
        
      }
      
      partial_derivs_vec[[k + counter + 2]] <- partial_deriv
      percent_models[k + counter + 2] <- 100*length(filtered_models)/numModels_filter_state
      
      filters_list[[k + counter + 2]] <- filtered_models
      
    } else {
      
      if ((partial_states_list[[k + counter + 2]][[k]][[1]] < current_state[[k]][[1]])) {
        
        bound <- 1
        
      } else {
        
        bound <- 2
        
      }
      
      partial_state_costs_object <- simulate_partial_state(partial_state = partial_states_list[[k + counter + 2]], current_state_rset_object = current_state_score_object,
                                                           reference_rset_object = reference_rset_object, partial_rset = partial_rset, scoring_method = scoring_method, hill_params = hill_params, param = k, bound = bound)
      
      if (bound == 1) {
        
        partial_deriv <- -(partial_state_costs_object$cf_vec_sum - current_state_score_object$cf_vec_sum)/update_vector_rep[k + counter + 2]
        
      } else {
        
        partial_deriv <- (partial_state_costs_object$cf_vec_sum - current_state_score_object$cf_vec_sum)/update_vector_rep[k + counter + 2]
        
      }
      
      partial_derivs_vec[k + counter + 2] <- partial_deriv
      percent_models[k + counter + 2] <- 100*nrow(partial_state_costs_object$partial_state_expression)/numModels_filter_state
      
    }
    ########################################
    if ((partial_states_list[[k + counter + 3]][[k]][[1]] >= current_state[[k]][[1]]) & (partial_states_list[[k + counter + 3]][[k]][[2]] <= current_state[[k]][[2]])) {
      
      #filter models of the current state using the param bounds of the partial state
      
      filtered_models <- which(current_state_params[,k] >= partial_states_list[[k + counter + 3]][[k]][[1]] & current_state_params[,k] <= partial_states_list[[k + counter + 3]][[k]][[2]])
      
      #round the partial states if the hill coefficient was changed
      if (k %in% hill_params) {
        
        filtered_models <- which(current_state_params[,k] >= round(partial_states_list[[k + counter + 3]][[k]][[1]]) & current_state_params[,k] <= round(partial_states_list[[k + counter + 3]][[k]][[2]]))
      }
      
      filtered_scores_object <- filter_score(current_state_rset_object = current_state_score_object, reference_rset_object = reference_rset_object, scoring_method = scoring_method, filter = filtered_models)
      
      if ((partial_states_list[[k + counter + 3]][[k]][[2]] < current_state[[k]][[2]])) {
        
        partial_deriv <- -(filtered_scores_object$cf_vec_sum - current_state_score_object$cf_vec_sum)/update_vector_rep[k + counter + 3]
        
      } else {
        
        partial_deriv <- (filtered_scores_object$cf_vec_sum - current_state_score_object$cf_vec_sum)/update_vector_rep[k + counter + 3]
        
      }
      
      partial_derivs_vec[[k + counter + 3]] <- partial_deriv
      percent_models[k + counter + 3] <- 100*length(filtered_models)/numModels_filter_state
      
      filters_list[[k + counter + 3]] <- filtered_models
      
    } else {
      
      if ((partial_states_list[[k + counter + 3]][[k]][[1]] < current_state[[k]][[1]])) {
        
        bound <- 1
        
      } else {
        
        bound <- 2
        
      }
      
      partial_state_costs_object <- simulate_partial_state(partial_state = partial_states_list[[k + counter + 3]], current_state_rset_object = current_state_score_object,
                                                           reference_rset_object = reference_rset_object, partial_rset = partial_rset, scoring_method = scoring_method, hill_params = hill_params, param = k, bound = bound)
      
      if (bound == 1) {
        
        partial_deriv <- -(partial_state_costs_object$cf_vec_sum - current_state_score_object$cf_vec_sum)/update_vector_rep[k + counter + 3]
        
      } else {
        
        partial_deriv <- (partial_state_costs_object$cf_vec_sum - current_state_score_object$cf_vec_sum)/update_vector_rep[k + counter + 3]
        
      }
      
      partial_derivs_vec[k + counter + 3] <- partial_deriv
      percent_models[k + counter + 3] <- 100*nrow(partial_state_costs_object$partial_state_expression)/numModels_filter_state
      
    }
    ########################################
    counter <- counter + 3
    
  }
  
  derivs_vec_minus <- vector(length = 2*number_of_params)
  derivs_vec_plus <- vector(length = 2*number_of_params)
  derivs_vec_mean <- vector(length = 2*number_of_params)
  
  counter <- 0
  for (z in 1:length(derivs_vec_minus)) {
    
    derivs_vec_minus[z] <- partial_derivs_vec[z + counter]
    
    counter <- counter + 1
    
  }
  
  counter <- 0
  z <- 1
  
  for (z in 1:length(derivs_vec_plus)) {
    
    derivs_vec_plus[z] <- partial_derivs_vec[z + counter + 1]
    
    counter <- counter + 1
    
  }
  
  z <- 1
  
  for (z in 1:length(derivs_vec_mean)) {
    
    derivs_vec_mean[z] <- mean(c(derivs_vec_minus[z], derivs_vec_plus[z]))
    
  }
  
  # taylor_derivs_vec <- vector(length = 2*number_of_params)
  # counter1 <- 0
  # counter2 <- 0
  # 
  # derivs_average <- vector(length = 2*number_of_params)
  # 
  # for (z in 1:number_of_params) {
  #   
  #   derivs_average[z] <- (partial_derivs_vec[z + counter2] + partial_derivs_vec[z + counter2 + 1])/(2*[z])
  #   derivs_average[z + 1] <- (partial_derivs_vec[z + counter2 + 2] + partial_derivs_vec[z + counter2 + 3])/(2*update_vector[z])
  #   
  #   taylor_derivs_vec[z + counter1] <- (partial_derivs_vec[z + counter2] + partial_derivs_vec[z + counter2 + 1] - 2*current_state_score_object$cf_vec_sum)/(2*update_vector[z])
  #   taylor_derivs_vec[z + 1 + counter1] <- (partial_derivs_vec[z + counter2 + 2] + partial_derivs_vec[z + counter2 + 3] - 2*current_state_score_object$cf_vec_sum)/(2*update_vector[z])
  #  
  #   counter1 <- counter1 + 1
  #   counter2 <- counter2 + 3
  #    
  # }
  # 
  
  derivs_object <- list(partial_derivs_vec = partial_derivs_vec, percent_models = percent_models, derivs_vec_plus = derivs_vec_plus, derivs_vec_minus = derivs_vec_minus, derivs_vec_mean = derivs_vec_mean)
  
  if (save_all == TRUE) {
    
    derivs_object <- list(partial_derivs_vec = partial_derivs_vec, percent_models = percent_models, derivs_vec_plus = derivs_vec_plus, derivs_vec_minus = derivs_vec_minus, derivs_vec_mean = derivs_vec_mean, partial_states_list = partial_states_list,
                          current_state_score_object = current_state_score_object, reference_rset_object = reference_rset_object, partial_rset = partial_rset, filters_list = filters_list, update_vector  = update_vector)
    
  }
  
  derivs_object
  
  
}
  ################################################################################################################

# #lfbgs, quasi newton, stochastic optimization
# a <- calc_partial_derivatives(current_state = ref_current_state0, reference_rset_object = reference_object, global_bounds = global_bounds_list, step_size = 0.05, numModels_filter_state = 2000, integrateStepSize_filter_state = 0.02, scoring_method = "Default",
#                                save_all = TRUE)
# a$partial_states_list[[5]]
# 
# b <- simulate_partial_state(partial_state = a$partial_states_list[[1]], current_state_rset_object = a$current_state_score_object, reference_rset_object = reference_object, partial_rset = a$partial_rset, scoring_method = "Default", hill_params = c(7,8), save_all = TRUE, param = 1, bound = 1)
# 
# test_partial_state <- a$partial_states_list[[1]]
# test_partial_state[[1]] <- list(1, 100)
# test_partial_state <- ref_state0
# 
# 
# c <- simulate_partial_state(partial_state = test_partial_state, current_state_rset_object = a$current_state_score_object, reference_rset_object = reference_object, partial_rset = a$partial_rset, scoring_method = "Default", hill_params = c(7,8), save_all = TRUE)
# 
# filtered_out <- 1:2000
# filtered_out <- filtered_out[-(a$filters_list[[2]])]
# 
# d <- filter_score(current_state_rset_object = a$current_state_score_object, reference_rset_object = reference_object, scoring_method = "Default", filter = filtered_out)
# e <- filter_score(current_state_rset_object = a$current_state_score_object, reference_rset_object = reference_object, scoring_method = "Default", filter = a$filters_list[[2]])
# 
# 
# plot(density(a$current_state_score_object$simulated_expression[, 1]))
# lines(density(d$filtered_expression[,1]), col = "red")
# lines(density(e$filtered_expression[,1]), col = "green")
# 
# 
# ###
# plot(density(sracipeParams(a$current_state_score_object$simulated_rset)[, 1]))
# min(sracipeParams(a$current_state_score_object$simulated_rset)[,1])
# 
# range(sracipeParams(a$current_state_score_object$simulated_rset)[,1])
# a$current_state_score_object$arguments$state
# 
# lines(density(sracipeParams(b$partial_state_expression_rset)[, 1]), col = "red")
# min(sracipeParams(b$partial_state_expression_rset)[, 1])
# 
# a$partial_states_list[[1]]
# 
# range(sracipeParams(a$current_state_score_object$simulated_rset)[,1])
# range(sracipeParams(b$partial_state_expression_rset)[, 1])
# 
# ###
# plot(density(a$current_state_score_object$simulated_expression[, 1]))
# lines(density(b$partial_state_expression[, 1]), col = "red")
# nrow(b$partial_state_expression)
# 
# ###
# min(sracipeParams(c$partial_state_expression_rset)[, 1])
# plot(density(a$current_state_score_object$simulated_expression[, 1]))
# lines(density(c$partial_state_expression[, 1]), col = "red")
# lines(density(c$partial_state_expression[, 1]), col = "red")
# nrow(c$partial_state_expression)
# c$cf_vec_sum
# ######
# 
# ###
# 
# 
# 
# 
# 
# 
# ###
# partial_state_expression_rset <- sracipeSimulate(partial_rset, integrate = TRUE, genParams = FALSE)
# 
# plot(density(log2(t(assay(a$current_state_score_object$simulated_rset)))))
# lines(density(log2(t(assay(c$partial_state_expression_rset)))), col = "red")
# 

#####################partial_derivs_loop##################################################################################################################################################################
loop_derivs <- function(number_of_iterations = 250, current_state, reference_rset_object, global_bounds, step_size = 0.01, numModels_filter_state = 20000, integrateStepSize_filter_state = 0.02, scoring_method = "Default") {
  
  partial_derivs_loop_list <- vector(mode = "list", length = number_of_iterations)
  
  for (i in 1:number_of_iterations) {
    
    partial_derivs_loop_list[[i]] <- calc_partial_derivatives(current_state = current_state, reference_rset_object = reference_rset_object, global_bounds = global_bounds, step_size = step_size,
                                                              numModels_filter_state = numModels_filter_state, integrateStepSize_filter_state = integrateStepSize_filter_state, scoring_method = scoring_method)
    
  }
                                 
  
  partial_derivs_loop_object <- list(partial_derivs_loop_list = partial_derivs_loop_list, arguments = list(current_state = current_state, reference_rset_object = reference_rset_object, global_bounds = global_bounds,
                                                                                                           step_size = step_size,
                                                                                                           numModels_filter_state = numModels_filter_state, integrateStepSize_filter_state = integrateStepSize_filter_state, scoring_method = scoring_method
                                                                                                           ))
  
}


#####################extract_statistics##################################################################################################################################################################
extract_statistics <- function(looped_derivs_object, reference_rset_object) {
  
  #vector of paramater names
  parameter_names <- colnames(sracipeParams(reference_rset_object$reference_rset))
  
  #vector of parameter_names + bound
  parameter_bounds_names <- paste(rep(parameter_names, each = 2),
                                  rep(c("lower", "upper"), times = length(parameter_names)), sep = "_")
  
  number_of_bounds <- 2*length(looped_derivs_object$arguments$current_state)
  number_of_iterations <- length(looped_derivs_object$partial_derivs_loop_list)
  
  
  partial_derivs_list_minus <- vector(mode = "list", length = number_of_bounds)
  partial_derivs_list_plus <- vector(mode = "list", length = number_of_bounds)
  partial_derivs_list_means <- vector(mode = "list", length = number_of_bounds)
  
  for (i in 1:number_of_bounds) {
    
    partials_vec_minus <- vector(length = number_of_iterations)
    partials_vec_plus <- vector(length = number_of_iterations)
    partials_vec_means <- vector(length = number_of_iterations)
    
    for (j in 1:number_of_iterations) {
      
      partials_vec_minus[j] <- looped_derivs_object$partial_derivs_loop_list[[j]]$derivs_vec_minus[i]
      partials_vec_plus[j] <- looped_derivs_object$partial_derivs_loop_list[[j]]$derivs_vec_plus[i]
      partials_vec_means[j] <- looped_derivs_object$partial_derivs_loop_list[[j]]$derivs_vec_mean[i]
      
    }
    
    partial_derivs_list_minus[[i]] <- partials_vec_minus
    partial_derivs_list_plus[[i]] <- partials_vec_plus
    partial_derivs_list_means[[i]] <- partials_vec_means
    
  }
  
  
  plots_list <-  vector(mode = "list", length = number_of_bounds)
  
  for (k in 1:number_of_bounds) {
    
    df <- data.frame(derivs_minus = partial_derivs_list_minus[[k]], derivs_plus = partial_derivs_list_plus[[k]], derivs_means = partial_derivs_list_means[[k]])
    
    plot <- ggplot(df, aes(x = derivs_minus)) + geom_density(color = "blue") + geom_density(aes(x = derivs_plus), color = "red")
    plot <- plot + geom_vline(aes(xintercept = mean(derivs_means)), color="black", linetype="dashed", size=1)
    plot <- plot + geom_vline(aes(xintercept = (mean(derivs_means) + sd(derivs_means))), color="purple", linetype="dashed", size=1) + geom_vline(aes(xintercept = (mean(derivs_means) - sd(derivs_means))), color="purple", linetype="dashed", size=1)
    plot <- plot + labs(title=parameter_bounds_names[k], x="Partial_Derivs", y = "Density")
  
    plots_list[[k]] <- plot
  
  }
  
  
  statistics_object <- list(partial_derivs_list_minus = partial_derivs_list_minus, partial_derivs_list_plus = partial_derivs_list_plus, partial_derivs_list_means = partial_derivs_list_means, plots_list = plots_list)
  statistics_object
  
}


# stats_test <- extract_statistics(looped_derivs_object = looped_derivs_01_10000_Default_current_state_between, reference_rset_object = reference_object)
# 
do.call("grid.arrange", stats_test$plots_list)
looped_derivs_01_10000_Default_current_state_between$arguments$current_state
looped_derivs_01_10000_Default_current_state_between$arguments$reference_rset_object$reference_state
#####################calc_partial_state##################################################################################################################################################################
partial_deriv_distribution <- function(partial_state, current_state, reference_rset_object, scoring_method, numModels = 10000, integrateStepSize = 0.02, param_number = 1, param_bound = 1, step_percent = 0.05, number_of_iterations, sim_current_state = TRUE) {
  
  partial_derivs_vec <- vector(length = number_of_iterations)
  
  current_state_score_object <- calc_score(state = current_state, reference_rset_object = reference_rset_object, scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, save_rset = TRUE)
  
  for (i in 1:number_of_iterations) {
    
    if (sim_current_state) {
  
      current_state_score_object <- calc_score(state = current_state, reference_rset_object = reference_rset_object, scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, save_rset = TRUE)
          
    }
  
  partial_state_score_object <- calc_score(state = partial_state, reference_rset_object = reference_rset_object, scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, save_rset = TRUE)
  
  step_size <- step_percent*(current_state[[param_number]][[2]] - current_state[[param_number]][[1]])
  
  
  if ((partial_state[[param_number]][[param_bound]] < current_state[[param_number]][[param_bound]])) {
    
    partial_deriv <- -(partial_state_score_object$cf_vec_sum - current_state_score_object$cf_vec_sum)/step_size
    
  } else {
    
    partial_deriv <- (partial_state_score_object$cf_vec_sum - current_state_score_object$cf_vec_sum)/step_size
    
  }
  
  partial_derivs_vec[i] <- partial_deriv
  
}
  
  
  
  
  partial_derivs_object <- list(partial_derivs_vec = partial_derivs_vec)
  
}
  
#####################filter_score_test##################################################################################################################################################################  

filter_score_test <- function(partial_state, current_state, reference_rset_object, param_number, param_bound, number_of_iterations = 200, numModels = 10000, integrateStepSize = 0.02, scoring_method = "Default", step_percent = 0.05) {
  
  
  current_state_score_object <- calc_score(state = current_state, reference_rset_object = reference_rset_object, scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, save_rset = TRUE)
  
  current_state_params <- sracipeParams(current_state_score_object$simulated_rset)
  
  
  filtered_models <- which(current_state_params[,param_number] >= partial_state[[param_number]][[1]] & current_state_params[,param_number] <= partial_state[[param_number]][[2]])
  # filtered_models2 <- intersect(which(current_state_params[,param_number] >= partial_state[[param_number]][[1]]), which(current_state_params[,param_number] <= partial_state[[param_number]][[2]]))
  # print(identical(filtered_models, filtered_models2))
  
  filtered_scores_object <- filter_score(current_state_rset_object = current_state_score_object, reference_rset_object = reference_rset_object, scoring_method = scoring_method, filter = filter)
  
  step_size <- step_percent*(current_state[[param_number]][[2]] - current_state[[param_number]][[1]])
  
  if ((partial_state[[param_number]][[param_bound]] < current_state[[param_number]][[param_bound]])) {
    
    partial_deriv <- -(filtered_scores_object$cf_vec_sum - current_state_score_object$cf_vec_sum)/step_size
    
  } else {
    
    partial_deriv <- (filtered_scores_object$cf_vec_sum - current_state_score_object$cf_vec_sum)/step_size
    
  }
  
  
}

# filter_score_test(partial_state = derivs_005_10000_Default_current_state_below$partial_states_list[[2]], current_state = current_state_below, reference_rset_object = reference_object2, scoring_method = "Default",
#                   param_number = 1, number_of_iterations = 200)


#RED
#Simulate current
#Filter models
#calc deriv


#Orange
#simulate current
#simulate partial
#calc deriv


#BLUE
#simulate current
#add models
#calc deriv

#PURPLE
#simulate
#simulate partial
#calc deriv

#GREEN
#simulaate current once
#simulate ppartial
#calc deriv



#####################deltaXplots##################################################################################################################################################################  
deltaXplot <- function(extracted_data_list, step_vector, reference_rset_object) {
    
    #vector of paramater names
    parameter_names <- colnames(sracipeParams(reference_rset_object$reference_rset))
    
    #vector of parameter_names + bound
    parameter_bounds_names <- paste(rep(parameter_names, each = 2),
                                    rep(c("lower", "upper"), times = length(parameter_names)), sep = "_")
    
  
  number_of_variables <- length(extracted_data_list[[1]]$partial_derivs_list_minus)
  
  score_vec_list <- vector(mode = "list", length = number_of_variables)
  
  for (i in 1:number_of_variables) {
    
    scores_vec <- vector(length = length(extracted_data_list))
   
    for (j in 1:length(extracted_data_list)) {
      
      ks_score <- dgof::ks.test(x = extracted_data_list[[j]]$partial_derivs_list_minus[[i]], y = extracted_data_list[[j]]$partial_derivs_list_plus[[i]])
      scores_vec[j] <- as.numeric(ks_score$statistic)
      
    }
    
    score_vec_list[[i]] <- scores_vec
     
  }
  
  plots_list <- vector(mode = "list", length = number_of_variables)
  
  for (k in 1:number_of_variables) {
    
    df <- data.frame(deltaX = step_vector, ks_score = score_vec_list[[k]])
    
    plot <- ggplot(data = df, aes(x = deltaX, y = ks_score)) + geom_line() + geom_point() + ylim(0, 1) + labs(title=parameter_bounds_names[k]) + geom_hline(aes(yintercept = 0.2), color = "red")
    
    plots_list[[k]] <- plot
    
    # print(df)
  }
  
  output <- list(score_vec_list = score_vec_list, plots_list = plots_list)
  output
  
}

########################################################################################################################################################################################################  


current_state_between <- list(list(20,50), list(20,50), list(0.2,0.5), list(0.2, 0.5), list(10,40), list(10,40), list(2,4), list(2,4), list(20,50), list(20,50))
ref_state0 <- list(list(1,100), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))

topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))

reference_object0 <- generate_ref_object(topology = topology_TS, reference_state = ref_state0, numModels = 10000, integrateStepSize = 0.02)

global_bounds_list <- list(G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0, TH_upper = 55, N_lower = 1, N_upper = 5, FC_lower = 1, FC_upper = 100)


looped_derivs_0001_10000_Default_current_state_between <- loop_derivs(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, global_bounds = global_bounds_list, step_size = 0.001, numModels_filter_state = 10000, integrateStepSize_filter_state = 0.02,
                                                                     scoring_method = "Default")

statistics_csbetween_ss0001 <- extract_statistics(looped_derivs_object = looped_derivs_0001_10000_Default_current_state_between, reference_rset_object = reference_object0)


looped_derivs_001_10000_Default_current_state_between <- loop_derivs(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, global_bounds = global_bounds_list, step_size = 0.01, numModels_filter_state = 10000, integrateStepSize_filter_state = 0.02,
                                                                     scoring_method = "Default")

statistics_csbetween_ss001 <- extract_statistics(looped_derivs_object = looped_derivs_001_10000_Default_current_state_between, reference_rset_object = reference_object0)

do.call("grid.arrange", statistics_csbetween_ss001$plots_list)


looped_derivs_005_10000_Default_current_state_between <- loop_derivs(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, global_bounds = global_bounds_list, step_size = 0.05, numModels_filter_state = 10000, integrateStepSize_filter_state = 0.02,
                                                                     scoring_method = "Default")

statistics_csbetween_ss005 <- extract_statistics(looped_derivs_object = looped_derivs_005_10000_Default_current_state_between, reference_rset_object = reference_object0)


looped_derivs_002_10000_Default_current_state_between <- loop_derivs(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, global_bounds = global_bounds_list, step_size = 0.02, numModels_filter_state = 10000, integrateStepSize_filter_state = 0.02,
                                                                    scoring_method = "Default")

looped_derivs_0005_10000_Default_current_state_between <- loop_derivs(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, global_bounds = global_bounds_list, step_size = 0.005, numModels_filter_state = 10000, integrateStepSize_filter_state = 0.02,
                                                                     scoring_method = "Default")

statistics_csbetween_ss0005 <- extract_statistics(looped_derivs_object = looped_derivs_0005_10000_Default_current_state_between, reference_rset_object = reference_object0)


looped_derivs_01_10000_Default_current_state_between <- loop_derivs(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, global_bounds = global_bounds_list, step_size = 0.1, numModels_filter_state = 10000, integrateStepSize_filter_state = 0.02,
                                                                      scoring_method = "Default")


stats_test0_0005 <- extract_statistics(looped_derivs_object = looped_derivs_0005_10000_Default_current_state_between, reference_rset_object = reference_object0)

stats_test0 <- extract_statistics(looped_derivs_object = looped_derivs_002_10000_Default_current_state_between, reference_rset_object = reference_object0)
do.call("grid.arrange", stats_test0$plots_list)

stats_test <- extract_statistics(looped_derivs_object = looped_derivs_01_10000_Default_current_state_between, reference_rset_object = reference_object0)
do.call("grid.arrange", stats_test$plots_list)

stats_test0_005 <- extract_statistics(looped_derivs_object = looped_derivs_005_10000_Default_current_state_between, reference_rset_object = reference_object0)
do.call("grid.arrange", stats_test0_005$plots_list)


deltaXscores_current_state_between <- deltaXplot(list(stats_test0_0005, stats_test0, stats_test0_005, stats_test), step_vector = c(0.005, 0.02, 0.05, 0.1), reference_rset_object = reference_object0)
do.call("grid.arrange", deltaXscores_current_state_between$plots_list)



###############
current_state_outside <- list(list(20,80), list(20,80), list(0.3,0.8), list(0.3, 0.8), list(10,45), list(10,45), list(2,6), list(2,6), list(20,80), list(20,80))
ref_state1 <- list(list(35,65), list(35,65), list(0.35,0.65), list(0.35, 0.65), list(15,40), list(15,40), list(3,5), list(3,5), list(35,65), list(35,65))


topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))
reference_object1 <- generate_ref_object(topology = topology_TS, reference_state = ref_state1, numModels = 10000, integrateStepSize = 0.02)
global_bounds_list <- list(G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0, TH_upper = 55, N_lower = 1, N_upper = 5, FC_lower = 1, FC_upper = 100)

looped_derivs_005_10000_Default_current_state_outside <- loop_derivs(number_of_iterations = 200, current_state = current_state_outside, reference_rset_object = reference_object1, global_bounds = global_bounds_list, step_size = 0.05, numModels_filter_state = 10000, integrateStepSize_filter_state = 0.02,
                                                                    scoring_method = "Default")

stats_test2 <- extract_statistics(looped_derivs_object = looped_derivs_005_10000_Default_current_state_outside, reference_rset_object = reference_object1)
do.call("grid.arrange", stats_test2$plots_list)

derivs_005_10000_Default_current_state_outside <- calc_partial_derivatives(current_state = current_state_outside, reference_rset_object = reference_object1, global_bounds = global_bounds_list, step_size = 0.05, numModels_filter_state = 10000, integrateStepSize_filter_state = 0.02, 
                         scoring_method = "Default", save_all = TRUE)

derivs_005_10000_Default_current_state_outside$partial_states_list[[1]]
##############

###############
current_state_below <- list(list(10,40), list(10,40), list(0.15,0.45), list(0.15, 0.45), list(5,20), list(5,20), list(2,3), list(2,3), list(10,40), list(10,40))

ref_state2 <- list(list(50,100), list(50,100), list(0.55,1), list(0.55, 1), list(30,55), list(30,55), list(4,5), list(4,5), list(50,100), list(50,100))


topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))
reference_object2 <- generate_ref_object(topology = topology_TS, reference_state = ref_state2, numModels = 10000, integrateStepSize = 0.02)
global_bounds_list <- list(G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0, TH_upper = 55, N_lower = 1, N_upper = 5, FC_lower = 1, FC_upper = 100)

looped_derivs_005_10000_Default_current_state_below <- loop_derivs(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, global_bounds = global_bounds_list, step_size = 0.05, numModels_filter_state = 10000, integrateStepSize_filter_state = 0.02,
                                                                     scoring_method = "Default")

stats_test3 <- extract_statistics(looped_derivs_object = looped_derivs_005_10000_Default_current_state_below, reference_rset_object = reference_object2)
do.call("grid.arrange", stats_test3$plots_list)

looped_derivs_002_10000_Default_current_state_below <- loop_derivs(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, global_bounds = global_bounds_list, step_size = 0.02, numModels_filter_state = 10000, integrateStepSize_filter_state = 0.02,
                                                                   scoring_method = "Default")

stats_test3_002 <- extract_statistics(looped_derivs_object = looped_derivs_002_10000_Default_current_state_below, reference_rset_object = reference_object2)
do.call("grid.arrange", stats_test3_002$plots_list)

looped_derivs_0005_10000_Default_current_state_below <- loop_derivs(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, global_bounds = global_bounds_list, step_size = 0.005, numModels_filter_state = 10000, integrateStepSize_filter_state = 0.02,
                                                                   scoring_method = "Default")

stats_test3_0005 <- extract_statistics(looped_derivs_object = looped_derivs_0005_10000_Default_current_state_below, reference_rset_object = reference_object2)


looped_derivs_0001_10000_Default_current_state_below <- loop_derivs(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, global_bounds = global_bounds_list, step_size = 0.001, numModels_filter_state = 10000, integrateStepSize_filter_state = 0.02,
                                                                    scoring_method = "Default")

stats_test3_0001 <- extract_statistics(looped_derivs_object = looped_derivs_0001_10000_Default_current_state_below, reference_rset_object = reference_object2)
do.call("grid.arrange", stats_test3_0005$plots_list)


looped_derivs_001_10000_Default_current_state_below <- loop_derivs(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, global_bounds = global_bounds_list, step_size = 0.01, numModels_filter_state = 10000, integrateStepSize_filter_state = 0.02,
                                                                    scoring_method = "Default")

stats_test3_001 <- extract_statistics(looped_derivs_object = looped_derivs_001_10000_Default_current_state_below, reference_rset_object = reference_object2)

deltaXscores_current_state_below <- deltaXplot(list(stats_test3_0005, stats_test3_002, stats_test3), step_vector = c(0.005, 0.02, 0.05))
do.call("grid.arrange", deltaXscores_current_state_below$plots_list)

derivs_005_10000_Default_current_state_below <- calc_partial_derivatives(current_state = current_state_below, reference_rset_object = reference_object2, global_bounds = global_bounds_list, step_size = 0.05, numModels_filter_state = 10000, integrateStepSize_filter_state = 0.02, 
                                                                           scoring_method = "Default", save_all = TRUE)

derivs_005_10000_Default_current_state_below$partial_states_list[[2]]


current_state_below_param1_bound1_subtract <- partial_deriv_distribution(partial_state = derivs_005_10000_Default_current_state_below$partial_states_list[[1]], current_state = current_state_below, reference_rset_object = reference_object2, scoring_method = "Default",
                           param_number = 1, param_bound = 1, step_percent = 0.05, number_of_iterations = 200)

current_state_below_param1_bound1_add <- partial_deriv_distribution(partial_state = derivs_005_10000_Default_current_state_below$partial_states_list[[2]], current_state = current_state_below, reference_rset_object = reference_object2, scoring_method = "Default",
                                                                         param_number = 1, param_bound = 1, step_percent = 0.05, number_of_iterations = 200)

filter_score_test(partial_state = derivs_005_10000_Default_current_state_below$partial_states_list[[1]], current_state = current_state_below, reference_rset_object = reference_object2, scoring_method = "Default",
                  param_number = 1, number_of_iterations = 200)


df <- data.frame(derivs_minus = current_state_below_param1_bound1_subtract$partial_derivs_vec, derivs_plus = -current_state_below_param1_bound1_add$partial_derivs_vec)

plot <- ggplot(df, aes(x = derivs_minus)) + geom_density(color = "blue") + geom_density(aes(x = derivs_plus), color = "red") 

plot <- stats_test3$plots_list[[1]]
plot <- plot + geom_density(data = df,aes(x = derivs_minus), color = "purple") + geom_density(data = df,aes(x = derivs_plus), color = "orange")
plot
##############


current_state_above <-  list(list(50,100), list(50,100), list(0.55,1), list(0.55, 1), list(30,55), list(30,55), list(4,5), list(4,5), list(50,100), list(50,100))
ref_state3 <- list(list(10,40), list(10,40), list(0.15,0.45), list(0.15, 0.45), list(5,20), list(5,20), list(2,3), list(2,3), list(10,40), list(10,40))



topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))
reference_object3 <- generate_ref_object(topology = topology_TS, reference_state = ref_state3, numModels = 10000, integrateStepSize = 0.02)
global_bounds_list <- list(G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0, TH_upper = 55, N_lower = 1, N_upper = 5, FC_lower = 1, FC_upper = 100)

looped_derivs_005_10000_Default_current_state_above <- loop_derivs(number_of_iterations = 200, current_state = current_state_above, reference_rset_object = reference_object3, global_bounds = global_bounds_list, step_size = 0.05, numModels_filter_state = 10000, integrateStepSize_filter_state = 0.02,
                                                                   scoring_method = "Default")

stats_test2 <- extract_statistics(looped_derivs_object = looped_derivs_005_10000_Default_current_state_outside, reference_rset_object = reference_object1)
do.call("grid.arrange", stats_test2$plots_list)
##############
# 
# 
# df <- data.frame(derivs_minus = stats_test$partial_derivs_list_minus[[1]], derivs_plus = stats_test$partial_derivs_list_plus[[1]], derivs_means = stats_test$partial_derivs_list_means[[1]])
# df_melted <- melt(df)
# 
# 
# plot <- ggplot(df, aes(x = derivs_minus)) + geom_density(color = "blue") + geom_density(aes(x = derivs_plus), color = "red")
# plot <- plot + geom_vline(aes(xintercept = mean(derivs_means)), color="black", linetype="dashed", size=1)
# plot <- plot + geom_vline(aes(xintercept = (mean(derivs_means) + sd(derivs_means))), color="purple", linetype="dashed", size=1) + geom_vline(aes(xintercept = (mean(derivs_means) - sd(derivs_means))), color="purple", linetype="dashed", size=1)
# plot#+ geom_vline(aes(xintercept = 0.005, color="black", linetype="dashed", size=1))
# 
# 
# plot
# 
# 
# taylor_derivs_list <- list(mode = "list", length = 20)
# 
# for (i in 1:20) {
#   
#   taylor_vec <- vector(length = 250)
#   
#   for (j in 1:250) {
#     
#     taylor_vec[j] <- looped_derivs_05_10000_Default$partial_derivs_loop_list[[j]]$taylor_derivs_vec[i]
#     
#     
#   }
#   taylor_derivs_list[[i]] <- taylor_vec
#   
#   
# }
# 
# #current: 20 - 40
# #ref state: 1- 100
# #green: 18,.45
# #red: 21.5, 40
# 
# df <- as.data.frame(matrix(partial_derivs_list[[1]], ncol = 1))
# colnames(df) <- "derivs"
# 
# p <- ggplot(df, aes(x = derivs)) + geom_density() + geom_vline(aes(xintercept= 0.005),color="blue", linetype="dashed", size=1) +geom_density() 
# p
# plot(density(partial_derivs_list[[3]]), col = "green",  xlim = c(-.01, .01))
# lines(density(partial_derivs_list[[4]]), col = "red")
# lines(density(partial_derivs_list[[2]]), col = "red")
# lines(density(partial_derivs_list[[2]]), col = "red")
# 
# 
# plot(density(partial_derivs_list[[39]]), col = "green",  xlim = c(-.005, .005))
# lines(density(partial_derivs_list[[40]]), col = "red")
# 
# for (i in 1:20) {
#   
#   plot()
#   
# }
# 
# plot(density(partial_derivs_list[[1]]), col = "green",  xlim = c(-.005, .005))
# lines(density(partial_derivs_list[[2]]), col = "red")
# 
# 
# plot(density(partial_derivs_list[[3]]), col = "green",  xlim = c(-.005, .005))
# lines(density(partial_derivs_list[[4]]), col = "red")
# 
# 
# plot(density(partial_derivs_list[[5]]), col = "green",  xlim = c(-.005, .005))
# lines(density(partial_derivs_list[[6]]), col = "red")
# 
# 
# plot(density(partial_derivs_list[[7]]), col = "green",  xlim = c(-.005, .005))
# lines(density(partial_derivs_list[[8]]), col = "red")
# 
# plot(density(partial_derivs_list[[9]]), col = "green",  xlim = c(-.5, .5))
# lines(density(partial_derivs_list[[10]]), col = "red")
# 
# plot(density(partial_derivs_list[[11]]), col = "green",  xlim = c(-.5, .5))
# lines(density(partial_derivs_list[[12]]), col = "red")
# 
# plot(density(partial_derivs_list[[7]]), col = "green",  xlim = c(-.005, .005))
# lines(density(partial_derivs_list[[8]]), col = "red")
# 
# plot(density(partial_derivs_list[[7]]), col = "green",  xlim = c(-.005, .005))
# lines(density(partial_derivs_list[[8]]), col = "red")
# 
