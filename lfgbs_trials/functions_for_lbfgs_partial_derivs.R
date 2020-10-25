#####################Globals##################################################################################################################################################################
current_state_og <- list(list(20,50), list(20,50), list(0.2,0.5), list(0.2, 0.5), list(10,40), list(10,40), list(2,4), list(2,4), list(20,50), list(20,50))
current_state_og <- split(unlist(current_state_og), ceiling(seq_along(unlist(current_state_og))/2))

current_state_og2 <- list(list(20,50), list(20,50), list(0.2,0.5), list(0.2, 0.5), list(10,40), list(10,40), list(1.95,4), list(2,4.05), list(20,50), list(20,50))


current_state <- list(list(20,50), list(20,50), list(0.2,0.5), list(0.2, 0.5), list(10,40), list(10,40), list(20,50), list(20,50))
state_vector <- unlist(current_state)
skip <- c("N")

reference_state <- list(list(1,100), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))
topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))
numModels <- 10000
integrateStepSize <- 0.02
reference <- generate_ref_object(topology = topology_TS, reference_state = reference_state, numModels = numModels, integrateStepSize = integrateStepSize)

global_bounds <- list(G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0, TH_upper = 55, N_lower = 1, N_upper = 5, FC_lower = 1, FC_upper = 100)
g_bounds <- list(G = c(1, 100), K = c(0.1, 1), TH = c(0, 55), N = c(1, 5), FC = c(1, 100))
scoring_method <- "Default"
step_size <- 0.02
gv_numModels_FC <- 2

iteration_pdf <- 0
iteration_cf <- 0

partial_derivs_input <- list()


gv_current_state_scores_object <- c()

#####################calculate_score##################################################################################################################################################################
calc_score <- function(state, reference_rset_object, scoring_method = "Default", numModels = 10000, integrateStepSize = 0.02, save_rset = FALSE) {
  
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
  if (scoring_method == "PCs") {
    
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
  
  costs_object <- list(cf_vec = cf_vec, cf_vec_sum = sum(cf_vec))
  
  if (save_rset == TRUE) {
    
    costs_object <- list(cf_vec = cf_vec, cf_vec_sum = sum(cf_vec), simulated_rset = simulated_rset,  simulated_expression = simulated_expression,
                         arguments = list(scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, state = state))
    
  }
  costs_object
}

# test3 <- calc_score(state = current_state_og2, reference_rset_object = reference, save_rset = TRUE)
# plot(density(sracipeParams(test3$simulated_rset)[,8]))
# 
# test4 <- calc_score(state = current_state_og, reference_rset_object = reference, save_rset = TRUE)
# plot(density(sracipeParams(test4$simulated_rset)[,8]))
#####################calculate_partial_derivatives##################################################################################################################################################################
calculate_partial_derivatives <- function(state, reference_rset_object, numModels = 10000, integrateStepSize = 0.02, numModels_FC = 1, scoring_method, global_bounds, step_size, skipped_params, current_state_scores_object) {
  
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

############################################################################################################################################################################################################################
# test <- calculate_partial_derivatives(state = current_state_og, reference_rset_object = reference, global_bounds = g_bounds, scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize,
#                                       step_size = step_size, numModels_FC = gv_numModels_FC)
# 
# 
# test2 <- partial_derivs_calc4(number_of_iterations = 1, current_state = current_state_og, reference_rset_object = reference, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 200, scoring_method = scoring_method, global_bounds = global_bounds,
#                               step_size = 0.02, repeat_current = FALSE)
# test$partial_states_list
# test2$partial_states_list
# test2$current_scores_list

#####################calculate_partial_derivatives_wrapper############################################################################################################################################################################################################################
calculate_partial_derivatives_wrapper <- function(vector) {
  
  print(vector)
  
  iteration_pdf <<- iteration_pdf + 1
  print(paste("Partial Derivative Function Iteration: ", iteration_pdf))
  
  state_list <- split(vector, ceiling(seq_along(vector)/2))
  
  partial_derivs_input <<- append(partial_derivs_input, list(state_list))
  
  partial_derivatives_object <- calculate_partial_derivatives(state = state_list, reference_rset_object = reference, global_bounds = g_bounds, scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize,
                                                       skipped_params = skip, step_size = step_size, numModels_FC = gv_numModels_FC, current_state_scores_object = gv_current_state_scores_object)
  
  output <- unname(partial_derivatives_object$partial_derivatives)
  output
  
}

calculate_partial_derivatives_wrapper(state_vector)
#####################loop_partial_derivatives##################################################################################################################################################################

loop_partial_derivatives <- function(state, number_of_iterations = 200) {
  
  partial_derivs_list <- vector(mode = "list", length = number_of_iterations)
  
  for (i in 1:number_of_iterations) {
    
    
    partial_derivatives_object <- calculate_partial_derivatives(state = state, reference_rset_object = reference, numModels = numModels, integrateStepSize = integrateStepSize, numModels_FC = numModels_FC, scoring_method = scoring_method, global_bounds = g_bounds, step_size = step_size)
    partial_derivs_list[[i]] <- partial_derivatives_object$partial_derivatives
  }
  
  partial_derivs_list
  
}

looped_derivs_cs_between_FC1_ss002 <- loop_partial_derivatives(state = current_state_og, number_of_iterations = 200)
#####################plot_statistics##################################################################################################################################################################

plot_statistics <- function(partial_derivs_list, plots_list) {
  
  partial_derivs_per_param <- vector(mode = "list", length = length(partial_derivs_list[[1]]))
  
  for (i in 1:length(partial_derivs_per_param)) {
    
    partial_derivs <- vector(length = length(partial_derivs_list))
    
    for (j in 1:length(partial_derivs)) {
      
      partial_derivs[[j]] <- partial_derivs_list[[j]][[i]]
      
    }
    
    partial_derivs_per_param[[i]] <- partial_derivs
    
  }
  
  plots <- vector(mode = "list", length = length(partial_derivs_per_param))
  
  for (i in 1:length(partial_derivs_per_param)) {
    
    df <- data.frame(derivs = partial_derivs_per_param[[i]])
    
    plot <- plots_list[[i]]
    plot <- plot + geom_density(data = df, aes(x = derivs), color = "green")
    
    plots[[i]] <- plot
  }
  
  plots
}

test <- plot_statistics(partial_derivs_list = looped_derivs_cs_between_FC1_ss002, plots_list = statistics_csbetween_ss002$plots_list)
test <- test[-c(13:16)]
do.call("grid.arrange", test)

#######################################################################################################################################################################################






# test <- calc_partial_derivatives(current_state = current_state_og, reference_rset_object = reference, global_bounds = global_bounds, step_size = 0.02,
#                                                      numModels_filter_state = 10000, integrateStepSize_filter_state = 0.02, scoring_method = scoring_method, save_all = T)
# 
# 
# test2 <- partial_derivs_calc4(number_of_iterations = 2, current_state = current_state_og, reference_rset_object = reference, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 200, scoring_method = scoring_method, global_bounds = global_bounds,
#                               step_size = 0.02, repeat_current = FALSE)
# test$partial_states_list
# test2$partial_derivs_list_all
# test2$current_scores_list
# calculate_partial_derivatives(state_vector = state_vector, reference_rset_object = reference, numModels_FC = 1, scoring_method = scoring_method, global_bounds = g_bounds, step_size = 0.02, skipped_params = skip)
# 
# partial_derivs_calc3 <- function(number_of_iterations, current_state, reference_rset_object, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 500, scoring_method, plots_list, global_bounds, step_size, repeat_current = FALSE) {
#   
#   if (step_size*numModels <= numModels_partial) {
#     
#     rep_cstate_number <- (numModels_partial/step_size)/numModels
#     
#   } else {
#     
#     numModels_partial <- step_size*numModels
#     rep_cstate_number <- 1
#     
#   }
#   
#   partial_derivs_object <- calc_partial_derivatives(current_state = current_state, reference_rset_object = reference_rset_object, global_bounds = global_bounds_list, step_size = step_size,
#                                                     numModels_filter_state = numModels, integrateStepSize_filter_state = integrateStepSize, scoring_method = scoring_method, save_all = T)
#   
#   update_vector <- rep(partial_derivs_object$update_vector, each = 2)
#   
#   ###
#   current_state_scores_object <- calc_score(state = current_state, reference_rset_object = reference_rset_object, scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, save_rset = TRUE)
#   current_state_expr_rep <- do.call(rbind, replicate(rep_cstate_number, current_state_scores_object$simulated_expression, simplify=FALSE))
#   ###
#   
#   a <- vector(length = 20)
#   counter1 <- 0
#   counter2 <- 0
#   
#   for (i in 1:10) {
#     
#     a[i + counter2] <- i + counter1
#     
#     a[i + 1 + counter2] <- i + counter1 + 3
#     
#     
#     counter1 <- counter1 + 3
#     counter2 <- counter2 + 1
#     
#   }
#   
#   partial_states_list <- partial_derivs_object$partial_states_list[a]
#   
#   
#   #################
#   i <- 1
#   counter1 <- 0
#   for (i in 1:10) {
#     
#     if (partial_states_list[[i + counter1]][[i]][[1]] < current_state[[i]][[1]]) {
#       
#       partial_states_list[[i + counter1]][[i]][[2]] <- current_state[[i]][[1]]
#       
#     }
#     
#     if (partial_states_list[[i + 1 + counter1]][[i]][[2]] > current_state[[i]][[2]]) {
#       
#       partial_states_list[[i + 1 + counter1]][[i]][[1]] <- current_state[[i]][[2]]
#       
#     }
#     
#     counter1 <- counter1 + 1
#     
#   }
#   
#   
#   
#   partial_derivs_list_L <- vector(mode = "list", length = number_of_iterations)
#   partial_derivs_list_U <- vector(mode = "list", length = number_of_iterations)
#   
#   i <- 1
#   for (i in 1:number_of_iterations) {
#     
#     if (repeat_current) {
#       
#       current_state_scores_object <- calc_score(state = current_state, reference_rset_object = reference_rset_object, scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, save_rset = TRUE)
#       current_state_expr_rep <- do.call(rbind, replicate(rep_cstate_number, current_state_scores_object$simulated_expression, simplify = FALSE))
#       
#     }
#     
#     
#     partial_derivs_vec_L <- c()
#     partial_derivs_vec_U <- c()
#     
#     for (j in 1:length(partial_states_list)) {
#       
#       partial_state_scores_object <- calc_score(state = partial_states_list[[j]], reference_rset_object = reference_rset_object, scoring_method = scoring_method, numModels = numModels_partial, integrateStepSize = integrateStepSize, save_rset = TRUE, update_vector)
#       
#       combined_expr <- rbind(current_state_expr_rep, partial_state_scores_object$simulated_expression)
#       
#       #####
#       
#       cf_vec <- vector(length = ncol(current_state_scores_object$simulated_expression))
#       
#       for (w in 1:length(cf_vec)) {
#         
#         ks_test <- dgof::ks.test(x = reference_rset_object$expression_ref[,w], y = combined_expr[,w])
#         cf_vec[w] <- as.numeric(ks_test$statistic)
#         cf_sum <-sum(cf_vec)
#       }
#       
#       
#       if (j%%2 != 0) {
#         
#         partial_deriv_minus<- -(cf_sum - current_state_scores_object$cf_vec_sum)/update_vector[j]
#         
#         partial_derivs_vec_L <- c(partial_derivs_vec_L, partial_deriv_minus)
#         
#         
#       } else {
#         
#         partial_deriv_plus <- (cf_sum - current_state_scores_object$cf_vec_sum)/update_vector[j]
#         
#         partial_derivs_vec_U <- c(partial_derivs_vec_U, partial_deriv_plus)
#         
#       }
#       
#     }
#     
#     partial_derivs_list_L[[i]] <- partial_derivs_vec_L
#     partial_derivs_list_U[[i]] <- partial_derivs_vec_U
#     
#   }
#   
#   ##################
#   partial_derivs_list_minus <- vector(mode = "list", length = 10)
#   partial_derivs_list_plus <- vector(mode = "list", length = 10)
#   
#   
#   for (k in 1:10) {
#     
#     vec_minus <- vector(length = number_of_iterations)
#     vec_plus <- vector(length = number_of_iterations)
#     
#     for (l in 1:number_of_iterations) {
#       
#       vec_minus[l] <- partial_derivs_list_L[[l]][[k]]
#       
#       vec_plus[l] <- partial_derivs_list_U[[l]][[k]]
#     }
#     
#     partial_derivs_list_minus[[k]] <- vec_minus
#     partial_derivs_list_plus[[k]] <- vec_plus
#     
#   }
#   
#   ##################
#   param <- rep(1:10,each = 2)
#   
#   k <- 1
#   
#   for (k in 1:length(plots_list)) {
#     
#     if (k %%2 != 0) {
#       
#       df <- data.frame(derivs = partial_derivs_list_minus[[param[k]]])
#       
#       plot <- plots_list[[k]]
#       
#       plot <- plot + geom_density(df, mapping = aes(x = derivs), color = "green")
#       
#       plots_list[[k]] <- plot
#       
#     } else {
#       
#       df <- data.frame(derivs = partial_derivs_list_plus[[param[k]]])
#       
#       plot <- plots_list[[k]]
#       
#       plot <- plot + geom_density(df, mapping = aes(x = derivs), color = "green")
#       
#       plots_list[[k]] <- plot
#       
#     }
#     
#   }
#   
#   partial_derivs_object <- list(plots_list = plots_list, partial_derivs_L_byiteration = partial_derivs_list_L, partial_derivs_U_byiteration = partial_derivs_list_U, partial_derivs_minus_byparam = partial_derivs_list_minus,
#                                 partial_derivs_plus_byparam = partial_derivs_list_plus, current_state_scores_object = current_state_scores_object, current_state_expr_rep = current_state_expr_rep, numModels_partial = numModels_partial,
#                                 rep_cstate_number = rep_cstate_number, partial_states_list = partial_states_list, update_vector = update_vector, combined_expr = combined_expr,
#                                 arguments = list(current_state = current_state, reference_rset_object = reference_rset_object, numModels = numModels, integrateStepSize = integrateStepSize, numModels_partial = numModels_partial,
#                                                  scoring_method = scoring_method, global_bounds = global_bounds_list, step_size = step_size, repeat_current = repeat_current))
#   
# }

partial_derivs_calc4 <- function(number_of_iterations, current_state, reference_rset_object, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 500, scoring_method, global_bounds, step_size, repeat_current = FALSE) {
  
  if (step_size*numModels <= numModels_partial) {
    
    rep_cstate_number <- (numModels_partial/step_size)/numModels
    
  } else {
    
    numModels_partial <- step_size*numModels
    rep_cstate_number <- 1
    
  }
  
  partial_derivs_object <- calc_partial_derivatives(current_state = current_state, reference_rset_object = reference_rset_object, global_bounds = global_bounds_list, step_size = step_size,
                                                    numModels_filter_state = numModels, integrateStepSize_filter_state = integrateStepSize, scoring_method = scoring_method, save_all = T)
  
  update_vector <- rep(partial_derivs_object$update_vector, each = 2)
  
  ###
  current_state_scores_object <- calc_score(state = current_state, reference_rset_object = reference_rset_object, scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, save_rset = TRUE)
  current_state_expr_rep <- do.call(rbind, replicate(rep_cstate_number, current_state_scores_object$simulated_expression, simplify=FALSE))
  ###
  
  a <- vector(length = 20)
  counter1 <- 0
  counter2 <- 0
  
  for (i in 1:10) {
    
    a[i + counter2] <- i + counter1
    
    a[i + 1 + counter2] <- i + counter1 + 3
    
    
    counter1 <- counter1 + 3
    counter2 <- counter2 + 1
    
  }
  
  partial_states_list <- partial_derivs_object$partial_states_list[a]
  
  
  #################
  i <- 1
  counter1 <- 0
  for (i in 1:10) {
    
    if (partial_states_list[[i + counter1]][[i]][[1]] < current_state[[i]][[1]]) {
      
      partial_states_list[[i + counter1]][[i]][[2]] <- current_state[[i]][[1]]
      
    }
    
    if (partial_states_list[[i + 1 + counter1]][[i]][[2]] > current_state[[i]][[2]]) {
      
      partial_states_list[[i + 1 + counter1]][[i]][[1]] <- current_state[[i]][[2]]
      
    }
    
    counter1 <- counter1 + 1
    
  }
  
  
  
  partial_derivs_list_L <- vector(mode = "list", length = number_of_iterations)
  partial_derivs_list_U <- vector(mode = "list", length = number_of_iterations)
  partial_derivs_list_all <- vector(mode = "list", length = number_of_iterations)
  partial_scores_list <- vector(mode = "list", length = number_of_iterations)
  current_scores_list <- vector(mode = "list", length = number_of_iterations)
  
  i <- 1
  for (i in 1:number_of_iterations) {
    
    if (repeat_current) {
      
      current_state_scores_object <- calc_score(state = current_state, reference_rset_object = reference_rset_object, scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, save_rset = TRUE)
      current_state_expr_rep <- do.call(rbind, replicate(rep_cstate_number, current_state_scores_object$simulated_expression, simplify = FALSE))
      
    }
    
    
    partial_derivs_vec_L <- c()
    partial_derivs_vec_U <- c()
    
    partial_derivatives <- vector(length = length(partial_states_list))
    partial_scores <- vector(length = length(partial_states_list))
    current_score <- vector(length = length(partial_states_list))
    
    for (j in 1:length(partial_states_list)) {
      
      partial_state_scores_object <- calc_score(state = partial_states_list[[j]], reference_rset_object = reference_rset_object, scoring_method = scoring_method, numModels = numModels_partial, integrateStepSize = integrateStepSize, save_rset = TRUE)
      
      combined_expr <- rbind(current_state_expr_rep, partial_state_scores_object$simulated_expression)
      
      #####
      
      cf_vec <- vector(length = ncol(current_state_scores_object$simulated_expression))
      
      for (w in 1:length(cf_vec)) {
        
        ks_test <- dgof::ks.test(x = reference_rset_object$expression_ref[,w], y = combined_expr[,w])
        cf_vec[w] <- as.numeric(ks_test$statistic)
        cf_sum <-sum(cf_vec)
      }
      
      partial_scores[j] <- cf_sum
      current_score[j] <- current_state_scores_object$cf_vec_sum
      
      if (j%%2 != 0) {
        
        partial_deriv_minus <- -(cf_sum - current_state_scores_object$cf_vec_sum)/update_vector[j]
        
        partial_derivs_vec_L <- c(partial_derivs_vec_L, partial_deriv_minus)
        
        partial_derivatives[j] <- partial_deriv_minus
        
      } else {
        
        partial_deriv_plus <- (cf_sum - current_state_scores_object$cf_vec_sum)/update_vector[j]
        
        partial_derivs_vec_U <- c(partial_derivs_vec_U, partial_deriv_plus)
        
        partial_derivatives[j] <- partial_deriv_plus
        
      }
      
    }
    
    partial_derivs_list_L[[i]] <- partial_derivs_vec_L
    partial_derivs_list_U[[i]] <- partial_derivs_vec_U
    
    partial_derivs_list_all[[i]] <- partial_derivatives
    partial_scores_list[[i]] <- partial_scores
    current_scores_list[[i]] <- current_score
    
    
  }
  
  ##################
  partial_derivs_list_minus <- vector(mode = "list", length = 10)
  partial_derivs_list_plus <- vector(mode = "list", length = 10)
  
  
  for (k in 1:10) {
    
    vec_minus <- vector(length = number_of_iterations)
    vec_plus <- vector(length = number_of_iterations)
    
    for (l in 1:number_of_iterations) {
      
      vec_minus[l] <- partial_derivs_list_L[[l]][[k]]
      
      vec_plus[l] <- partial_derivs_list_U[[l]][[k]]
    }
    
    partial_derivs_list_minus[[k]] <- vec_minus
    partial_derivs_list_plus[[k]] <- vec_plus
    
  }
  
 
  partial_derivs_object <- list(update_vector = update_vector, partial_derivs_list_all = partial_derivs_list_all, partial_scores_list = partial_scores_list, current_scores_list = current_scores_list, partial_states_list = partial_states_list,
                                arguments = list(current_state = current_state, reference_rset_object = reference_rset_object, numModels = numModels, integrateStepSize = integrateStepSize, numModels_partial = numModels_partial,
                                                 scoring_method = scoring_method, global_bounds = global_bounds, step_size = step_size, repeat_current = repeat_current))
  
}




