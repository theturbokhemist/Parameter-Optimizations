eq1 <- function(x) {
  
  output <- 0.5*1*x
  output
  
}

curve(expr = eq2, from = 100, to = 1000)

eq2 <- function(x) {
  
  output <- 0.5*1*x/(0.5*x+1)
  output
  
}

curve(expr = eq1, from = 1, to = 100)
####################################Generate Update States########################################################################################################################################################################################################################
generate_update_states <- function(state_current, reference_rset_object, global_bounds, finite_difference = 0.1, skipped_params) {
  
  #######Pre-Processing######
  ##########################################transform state_reference into list

  state_reference <- split(unlist(reference_rset_object$reference_state), ceiling(seq_along(unlist(reference_rset_object$reference_state))/2))
  
  ##########################################extract parameter names and identify index of each parameter type
  parameter_names <- colnames(sracipeParams(reference_rset_object$reference_rset))
  
  parameter_bounds_names <- paste(rep(parameter_names, each = 2),
                                  rep(c("lower", "upper"), times = length(parameter_names)), sep = "_")
  
  production_params <- grep(pattern = "G", parameter_names)
  
  degradation_params <- grep(pattern = "K", parameter_names)
  
  threshold_params <- grep(pattern = "TH", parameter_names)
  
  hill_params <- grep(pattern = "N", parameter_names)
  
  FC_params <- grep(pattern = "FC", parameter_names)
  
  params_list <- list(production_params = production_params, degradation_params = degradation_params, threshold_params = threshold_params, hill_params = hill_params, FC_params = FC_params)
  
}

current_state_og <- list(c(20,50), c(20,50), c(0.2,0.5), c(0.2, 0.5), c(10,40), c(10,40), c(2,4), c(2,4), c(20,50), c(20,50))
current_state_og
generate_update_states(state_current = current_state_og, reference_rset_object = reference, global_bounds = g_bounds, finite_difference = 0.1, )
####################################Broydens Method########################################################################################################################################################################################################################

Broyden <- function(state, reference_rset_object, numModels = 10000, integrateStepSize = 0.02, numModels_FC = 1, scoring_method, global_bounds, step_size, skipped_params, finite_difference = 0.1) {
  
  #######Pre-Processing######
  ##########################################transform state_vector into list
  state_current <- state
  
  state_reference <- split(unlist(reference_rset_object$reference_state), ceiling(seq_along(unlist(reference_rset_object$reference_state))/2))
  
  ##########################################extract parameter names and identify index of each parameter type
  parameter_names <- colnames(sracipeParams(reference_rset_object$reference_rset))
  
  parameter_bounds_names <- paste(rep(parameter_names, each = 2),
                                  rep(c("lower", "upper"), times = length(parameter_names)), sep = "_")
  
  production_params <- grep(pattern = "G", parameter_names)
  
  degradation_params <- grep(pattern = "K", parameter_names)
  
  threshold_params <- grep(pattern = "TH", parameter_names)
  
  hill_params <- grep(pattern = "N", parameter_names)
  
  FC_params <- grep(pattern = "FC", parameter_names)
  
  params_list <- list(production_params = production_params, degradation_params = degradation_params, threshold_params = threshold_params, hill_params = hill_params, FC_params = FC_params)
  
  ##########################################index of parameters which should be skipped from full parameter list
  if (!missing(skipped_params)) {
    
    params_skipped <- c()
    bounds_skipped <- c()
    
    for (i in 1:length(skipped_params)) {
      
      params_skipped <- c(params_skipped, grep(pattern = skipped_params[i], parameter_names) )
      
    }
    
    for (i in 1:length(skipped_params)) {
      
      bounds_skipped <- c(bounds_skipped, grep(pattern = skipped_params[i], parameter_bounds_names) )
      
    }
    
  }
  
  ##########################################append reference values for skipped parameters
  if (!missing(skipped_params)) {
    
    if (sum(skipped_params %in% c("G")) > 0) {
      
      state_current <- append(state_current, state_reference[production_params], after = production_params[1] - 1)   
      
    }
    
    if (sum(skipped_params %in% c("K")) > 0) {
      
      state_current <- append(state_current, state_reference[degradation_params], after = degradation_params[1] - 1)   
      
    }
    
    if (sum(skipped_params %in% c("TH")) > 0) {
      
      state_current <- append(state_current, state_reference[threshold_params], after = threshold_params[1] - 1)   
      
    }
    
    if (sum(skipped_params %in% c("N")) > 0) {
      
      state_current <- append(state_current, state_reference[hill_params], after = hill_params[1] - 1)   
      
    }
    
    if (sum(skipped_params %in% c("FC")) > 0) {
      
      state_current <- append(state_current, state_reference[FC_params], after = FC_params[1] - 1)   
      
    }
  }
  
  names(state_current) <- parameter_names
  
  ##########################################global bounds for each parameter
  global_bounds_per_parameter <- list()
  
  for (i in 1:length(global_bounds)) {
    
    global_bounds_per_parameter <- append(global_bounds_per_parameter, (rep(global_bounds[i], length(params_list[[i]]))))
    
  }
  
  names(global_bounds_per_parameter) <- parameter_names
  
  # if (!missing(skipped_params)) {
  #   
  #   global_bounds_per_parameter <- global_bounds_per_parameter[-params_skipped]
  #   
  # }
  
  
  ##########################################check if current state bounds are within global bounds for each paramter
  
  for (i in 1:length(state_current)) {
    
    if (state_current[[i]][[1]] < global_bounds_per_parameter[[i]][[1]]) {
      
      state_current[[i]][[1]] <- global_bounds_per_parameter[[i]][[1]]
      
    }
    
    if (state_current[[i]][[2]] > global_bounds_per_parameter[[i]][[2]]) {
      
      state_current[[i]][[2]] <- global_bounds_per_parameter[[i]][[2]]
      
    }
    
  }
  
  #######Generate partial states list######
  
  ##########################################number of models to simulate to add
  numModels_partial <- numModels*step_size*numModels_FC
  
  ##########################################generate vector of the step magnitude for each parameter
  step_magnitudes <- vector(length = length(parameter_names))
  
  counter <- 0
  
  for (i in 1:length(step_magnitudes)) {
    
    step_magnitudes[i] <- (state_current[[i]][[2]] - state_current[[i]][[1]])*step_size
    
  }
  
  step_magnitudes_rep <- rep(step_magnitudes, each = 2)
  
  ##########################################generate list of partial states
  partial_states_list <- vector(mode = "list", length = length(parameter_names)*2) #*2 because subtracting from lower bound and adding to upper bound
  
  counter <- 0
  
  for (j in 1:length(state_current)) {
    
    #first lower/upper refers to which bound, second minus/plus refers t
    
    partial_state_lower_bound <- state_current
    partial_state_upper_bound <- state_current
    
    partial_state_lower_bound[[j]][[1]] <- state_current[[j]][[1]] - step_magnitudes[j]
    partial_state_lower_bound[[j]][[2]] <- state_current[[j]][[1]]
    
    partial_state_upper_bound[[j]][[1]] <- state_current[[j]][[2]]
    partial_state_upper_bound[[j]][[2]] <- state_current[[j]][[2]] + step_magnitudes[j]
    
    
    partial_states_list[[j + counter]] <- partial_state_lower_bound
    partial_states_list[[j + counter + 1]] <- partial_state_upper_bound
    
    counter <- counter + 1
    
  }
}
  ############################################################################################################################################################################################################################################################


hessian <- function(state, reference_rset_object, numModels = 10000, integrateStepSize = 0.02, numModels_FC = 1, scoring_method, global_bounds, step_size, skipped_params, current_state_scores_object) {
  
  #######Pre-Processing######
  ##########################################transform state_vector into list
  state_current <- state
  
  state_reference <- split(unlist(reference_rset_object$reference_state), ceiling(seq_along(unlist(reference_rset_object$reference_state))/2))
  
  ##########################################extract parameter names and identify index of each parameter type
  parameter_names <- colnames(sracipeParams(reference_rset_object$reference_rset))
  
  parameter_bounds_names <- paste(rep(parameter_names, each = 2),
                                  rep(c("lower", "upper"), times = length(parameter_names)), sep = "_")
  
  production_params <- grep(pattern = "G", parameter_names)
  
  degradation_params <- grep(pattern = "K", parameter_names)
  
  threshold_params <- grep(pattern = "TH", parameter_names)
  
  hill_params <- grep(pattern = "N", parameter_names)
  
  FC_params <- grep(pattern = "FC", parameter_names)
  
  params_list <- list(production_params = production_params, degradation_params = degradation_params, threshold_params = threshold_params, hill_params = hill_params, FC_params = FC_params)
  
  ##########################################index of parameters which should be skipped from full parameter list
  if (!missing(skipped_params)) {
    
    params_skipped <- c()
    bounds_skipped <- c()
    
    for (i in 1:length(skipped_params)) {
      
      params_skipped <- c(params_skipped, grep(pattern = skipped_params[i], parameter_names) )
      
    }
    
    for (i in 1:length(skipped_params)) {
      
      bounds_skipped <- c(bounds_skipped, grep(pattern = skipped_params[i], parameter_bounds_names) )
      
    }
    
  }
  
  ##########################################append reference values for skipped parameters
  if (!missing(skipped_params)) {
    
    if (sum(skipped_params %in% c("G")) > 0) {
      
      state_current <- append(state_current, state_reference[production_params], after = production_params[1] - 1)   
      
    }
    
    if (sum(skipped_params %in% c("K")) > 0) {
      
      state_current <- append(state_current, state_reference[degradation_params], after = degradation_params[1] - 1)   
      
    }
    
    if (sum(skipped_params %in% c("TH")) > 0) {
      
      state_current <- append(state_current, state_reference[threshold_params], after = threshold_params[1] - 1)   
      
    }
    
    if (sum(skipped_params %in% c("N")) > 0) {
      
      state_current <- append(state_current, state_reference[hill_params], after = hill_params[1] - 1)   
      
    }
    
    if (sum(skipped_params %in% c("FC")) > 0) {
      
      state_current <- append(state_current, state_reference[FC_params], after = FC_params[1] - 1)   
      
    }
  }
  
  names(state_current) <- parameter_names
  
  ##########################################global bounds for each parameter
  global_bounds_per_parameter <- list()
  
  for (i in 1:length(global_bounds)) {
    
    global_bounds_per_parameter <- append(global_bounds_per_parameter, (rep(global_bounds[i], length(params_list[[i]]))))
    
  }
  
  names(global_bounds_per_parameter) <- parameter_names
  
  # if (!missing(skipped_params)) {
  #   
  #   global_bounds_per_parameter <- global_bounds_per_parameter[-params_skipped]
  #   
  # }
  
  
  ##########################################check if current state bounds are within global bounds for each paramter
  
  for (i in 1:length(state_current)) {
    
    if (state_current[[i]][[1]] < global_bounds_per_parameter[[i]][[1]]) {
      
      state_current[[i]][[1]] <- global_bounds_per_parameter[[i]][[1]]
      
    }
    
    if (state_current[[i]][[2]] > global_bounds_per_parameter[[i]][[2]]) {
      
      state_current[[i]][[2]] <- global_bounds_per_parameter[[i]][[2]]
      
    }
    
  }
  
  #######Generate partial states list######
  
  ##########################################number of models to simulate to add
  numModels_partial <- numModels*step_size*numModels_FC
  
  ##########################################generate vector of the step magnitude for each parameter
  step_magnitudes <- vector(length = length(parameter_names))
  
  counter <- 0
  
  for (i in 1:length(step_magnitudes)) {
    
    step_magnitudes[i] <- (state_current[[i]][[2]] - state_current[[i]][[1]])*step_size
    
  }
  
  step_magnitudes_rep <- rep(step_magnitudes, each = 2)
  
  ##########################################generate list of partial states
  partial_states_list <- vector(mode = "list", length = length(parameter_names)*2) #*2 because subtracting from lower bound and adding to upper bound
  
  counter <- 0
  
  for (j in 1:length(state_current)) {
    
    #first lower/upper refers to which bound, second minus/plus refers t
    
    partial_state_lower_bound <- state_current
    partial_state_upper_bound <- state_current
    
    partial_state_lower_bound[[j]][[1]] <- state_current[[j]][[1]] - step_magnitudes[j]
    partial_state_lower_bound[[j]][[2]] <- state_current[[j]][[1]]
    
    partial_state_upper_bound[[j]][[1]] <- state_current[[j]][[2]]
    partial_state_upper_bound[[j]][[2]] <- state_current[[j]][[2]] + step_magnitudes[j]
    
    
    partial_states_list[[j + counter]] <- partial_state_lower_bound
    partial_states_list[[j + counter + 1]] <- partial_state_upper_bound
    
    counter <- counter + 1
    
  }
  
  #######Simulate partial states and Calculate Partial Derivatives######
  ##########################################Simulate current state
  
  # current_state_scores_object <- calc_score(state = state_current, reference_rset_object = reference_rset_object, scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, save_rset = TRUE)
  
  ##########################################Duplicate current state expression matrix *FC* times
  current_state_expr_rep <- do.call(rbind, replicate(numModels_FC, current_state_scores_object$simulated_expression, simplify = FALSE))
  
  ##########################################initialize partial derivatives vector
  partial_derivatives <- vector(length = length(partial_states_list))
  partial_scores <- vector(length = length(partial_states_list))
  current_score <- vector(length = length(partial_states_list))
  
  # partial_derivs_vec_L <- c()
  # partial_derivs_vec_U <- c()
  
  ##########################################loop through partial states list
  for (j in 1:length(partial_states_list)) {
    
    ##########################################Simulate partial states
    partial_state_scores_object <- calc_score(state = partial_states_list[[j]], reference_rset_object = reference_rset_object, scoring_method = scoring_method, numModels = numModels_partial, integrateStepSize = integrateStepSize, save_rset = TRUE)
    
    ##########################################Combine expression of simulated partial state with expression of the duplicated current state expression
    combined_expr <- rbind(current_state_expr_rep, partial_state_scores_object$simulated_expression)
    
    ##########################################Calculate ks score of the above combined expression
    cf_vec <- vector(length = ncol(combined_expr))
    
    for (w in 1:length(cf_vec)) {
      
      ks_test <- dgof::ks.test(x = reference_rset_object$expression_ref[,w], y = combined_expr[,w])
      cf_vec[w] <- as.numeric(ks_test$statistic)
      cf_sum <-sum(cf_vec)
    }
    
    partial_scores[j] <- cf_sum
    current_score[j] <- current_state_scores_object$cf_vec_sum
    ##########################################Calculate partial derivatives
    if (j%%2 != 0) {
      
      partial_derivatives[j] <- -(cf_sum - current_state_scores_object$cf_vec_sum)/step_magnitudes_rep[j]
      
      # partial_deriv_minus<- -(cf_sum - current_state_scores_object$cf_vec_sum)/update_vector[j]
      # 
      # partial_derivs_vec_L <- c(partial_derivs_vec_L, partial_deriv_minus)
      
      
    } else {
      
      partial_derivatives[j] <- (cf_sum - current_state_scores_object$cf_vec_sum)/step_magnitudes_rep[j]
      
      # partial_deriv_plus <- (cf_sum - current_state_scores_object$cf_vec_sum)/update_vector[j]
      # 
      # partial_derivs_vec_U <- c(partial_derivs_vec_U, partial_deriv_plus)
      
    }
    
  }
  
  names(partial_derivatives) <- parameter_bounds_names
  
  
  ##########################################remove "skipped" parameters
  
  if (!missing(skipped_params)) {
    
    partial_derivatives <- partial_derivatives[-bounds_skipped]
    
  }
  
  output_object <- list(partial_derivatives = partial_derivatives, partial_scores = partial_scores, current_score = current_score, sign = rep(c(-1, 1), length(parameter_names)), step_magnitudes = step_magnitudes_rep, partial_states_list = partial_states_list)
  output_object
  
}
