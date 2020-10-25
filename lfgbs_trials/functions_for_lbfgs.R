#####################Globals##################################################################################################################################################################
# current_state <- list(list(20,50), list(20,50), list(0.2,0.5), list(0.2, 0.5), list(10,40), list(10,40), list(2,4), list(2,4), list(20,50), list(20,50))
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

iteration_cf <- 0
cost_function_input <- list()

#####################cost function wrapper############################################################################################################################################################################################################################
cost_function_wrapper <- function(vector) {
  
  print(vector)
  
  iteration_cf <<- iteration_cf + 1
  print(paste("Cost Function Iteration: ", iteration_cf))
  
  state_list <- split(vector, ceiling(seq_along(vector)/2))
  
  cost_function_input <<- append(cost_function_input, list(state_list))
  
  costs_object <- cost_function(state = state_list, reference_rset_object = reference, global_bounds = g_bounds, scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, save_rset = FALSE, skipped_params = skip)
  
  output <- costs_object$cf_vec_sum
  
  output
  
}

#####################cost_function############################################################################################################################################################################################################################
cost_function <- function(state, reference_rset_object, global_bounds, scoring_method = "Default", numModels = 10000, integrateStepSize = 0.02, save_rset = FALSE, skipped_params) {
  
  #######Pre-Processing######
  ##########################################transform state_vector into list
  state_current <- state
  
  state_reference <- split(unlist(reference_rset_object$reference_state), ceiling(seq_along(unlist(reference_rset_object$reference_state))/2))
  
  ##########################################extract parameter names and identify index of each parameter type
  parameter_names <- colnames(sracipeParams(reference_rset_object$reference_rset))
  
  production_params <- grep(pattern = "G", parameter_names)
  
  degradation_params <- grep(pattern = "K", parameter_names)
  
  threshold_params <- grep(pattern = "TH", parameter_names)
  
  hill_params <- grep(pattern = "N", parameter_names)
  
  FC_params <- grep(pattern = "FC", parameter_names)
  
  params_list <- list(production_params = production_params, degradation_params = degradation_params, threshold_params = threshold_params, hill_params = hill_params, FC_params = FC_params)

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
  
  ##########################################check if current state bounds are within global bounds for each paramter
  
  for (i in 1:length(state_current)) {
    
    if (state_current[[i]][[1]] < global_bounds_per_parameter[[i]][[1]]) {
      
      state_current[[i]][[1]] <- global_bounds_per_parameter[[i]][[1]]
      
    }
    
    if (state_current[[i]][[2]] > global_bounds_per_parameter[[i]][[2]]) {
      
      state_current[[i]][[2]] <- global_bounds_per_parameter[[i]][[2]]
      
    }
    
  }
   
  
  #######Pre-Simulation and Simulation######
  ##########################################create rset object
  state <- state_current
  
  simulated_rset <- sracipeSimulate(reference_rset_object$topology, integrate = FALSE,
                                    numModels = numModels, genParams = TRUE,
                                    integrateStepSize = integrateStepSize)
  
  ##########################################sample parameter ranges from current state *numModels* times, add those results into the rset objects parameter values matrix
  for (t in 1:length(state)) {
    
    sracipeParams(simulated_rset)[,t] <- runif(numModels, min = state[[t]][[1]], max = state[[t]][[2]])  
    
  }
  
  #sample and round all the hill coefficients
  for (u in 1:length(hill_params)) {
    #you need to subtract 0.5 from lower bound and add 0.5 to upper bound in order to get a uniform distribution
    sracipeParams(simulated_rset)[,hill_params[u]] <- round(runif(numModels, min= round(state[[hill_params[u]]][[1]]) - 0.5, max= round(state[[hill_params[u]]][[2]]) + 0.5))
    
  }
  
  ##########################################simulate
  simulated_rset <- sracipeSimulate(simulated_rset, integrate = TRUE, genParams = FALSE)
  
  #######Calculations######
  
  ##########################################log transform
  simulated_rset_expr_log <- log2(t(assay(simulated_rset)))
  
  ##########################################z-score normalize
  simulated_rset_expr_log_z <- sweep(simulated_rset_expr_log,
                                     2, reference_rset_object$means_ref, FUN = "-")
  
  simulated_rset_expr_log_z <- sweep(simulated_rset_expr_log_z,
                                     2, reference_rset_object$sds_ref, FUN = "/")
  
  ##########################################rotate (transformation)
  simulated_rset_expr_PCs <- simulated_rset_expr_log_z %*% reference_rset_object$prcomp_ref$rotation
  
  #######Scoring######
  
  #cost_function vector
  cf_vec <- vector(length = ncol(simulated_rset_expr_log_z))
  
  ##########################################kolmogorov-smirnov
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
    
    costs_object <- list(cf_vec = cf_vec, cf_vec_sum = sum(cf_vec), simulated_rset = simulated_rset, simulated_expression = simulated_expression,
                         arguments = list(scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, state = state, params_list = params_list, reference_state = reference_state, reference_rset_object = reference_rset_object))
    
  }
  
  gv_current_state_scores_object <<- list(cf_vec = cf_vec, cf_vec_sum = sum(cf_vec), simulated_rset = simulated_rset, simulated_expression = simulated_expression,
                                          arguments = list(scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, state = state, params_list = params_list, reference_state = reference_state, reference_rset_object = reference_rset_object))
  
  costs_object
  
  
}


cost_function_wrapper(vector = state_vector)

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
  