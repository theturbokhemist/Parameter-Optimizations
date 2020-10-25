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

# calc_score(state = ref_state0, reference_rset_object = test, scoring_method = "PCs", numModels = 10000, integrateStepSize = 0.02)
##################################################################################################################################################################

ref_state0 <- list(list(1,100), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))
ref_current_state0 <- list(list(20,50), list(20,50), list(0.2,0.5), list(0.2, 0.5), list(10,40), list(10,40), list(2,4), list(2,4), list(20,50), list(20,50))
ref_current_state1 <- list(list(1,100), list(20,40), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))

topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))
scoring_method1 <- "Default"
reference_object <- generate_ref_object(topology = topology_TS, reference_state = ref_state0, numModels = 10000, integrateStepSize = 0.02)
global_bounds_list <- list(G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0, TH_upper = 55, N_lower = 1, N_upper = 5, FC_lower = 1, FC_upper = 100)

#####################calc_partial_derivs##################################################################################################################################################################
calc_partial_derivs <- function(current_state, reference_rset_object, global_bounds, step_size = 0.01, numModels = 10000, integrateStepSize = 0.02, numModels_filter_state = 20000, integrateStepSize_filter_state = 0.02, scoring_method = scoring_method,
                                use_global_bounds = FALSE, full_save = FALSE, current_state_costs_object) {
  
  filter_state <- current_state
  names(filter_state) <- colnames(sracipeParams(reference_rset_object$reference_rset))
  
  #number of parameters
  number_of_params <- length(current_state)
  
  #vector of paramater names
  parameter_names <-names(filter_state)
  
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
  
  if (use_global_bounds == TRUE) {
    
    for (i in 1:number_of_params) {
      
      if (i %in% production_params) {
        
        update_vector[i] <- (global_bounds$G_upper - global_bounds$G_lower)*step_size
        
        update_vector_rep[(i + counter):(i+counter+3)] <- rep(update_vector[i], 4)
      }
      
      if (i %in% degradation_params) {
        
        update_vector[i] <- (global_bounds$K_upper - global_bounds$K_lower)*step_size
        
        update_vector_rep[(i + counter):(i+counter+3)] <- rep(update_vector[i], 4)
      }
      
      if (i %in% threshold_params) {
        
        update_vector[i] <- (global_bounds$TH_upper - global_bounds$TH_lower)*step_size
        
        update_vector_rep[(i + counter):(i+counter+3)] <- rep(update_vector[i], 4)
      }
      
      if (i %in% hill_params) {
        
        update_vector[i] <- (global_bounds$N_upper - global_bounds$N_lower)*step_size
        
        update_vector_rep[(i + counter):(i+counter+3)] <- rep(update_vector[i], 4)
      }
      
      if (i %in% FC_params) {
        
        update_vector[i] <- (global_bounds$FC_upper - global_bounds$FC_lower)*step_size
        
        update_vector_rep[(i + counter):(i+counter+3)] <- rep(update_vector[i], 4)
      }
      counter <- counter + 3
    }
    
  
    } else {
    
      for (i in 1:number_of_params) {
          
          update_vector[i] <- (current_state[[i]][[2]] - current_state[[i]][[1]])*step_size
          
          update_vector_rep[(i + counter):(i+counter+3)] <- rep(update_vector[i], 4)
          
          counter <- counter + 3
          
        }
  }

  
  #generate list of partial states
  partial_states_list <- vector(mode = "list", length = number_of_params*4) #*2 for 2 bounds and *2 for adding and subtracting from each bound
  
  counter <- 0
  
  for (j in 1:length(filter_state)) {
    
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
  
  #generate filter state
  for (k in 1:length(filter_state)) {
    
    filter_state[[k]][1] <- current_state[[k]][[1]] - update_vector[k]
    filter_state[[k]][2] <- current_state[[k]][[2]] + update_vector[k]
    
  }
  
  #simulate filter state
  simulated_filter_state <- calc_score(state = filter_state, reference_rset_object = reference_rset_object, scoring_method = scoring_method, numModels = numModels_filter_state, integrateStepSize = integrateStepSize_filter_state, save_rset = TRUE)
  simulated_filter_state_params <- sracipeParams(simulated_filter_state$simulated_rset)

  
  #create list of filters
  partial_states_filters_list <- vector(mode = "list", length = length(partial_states_list))
  filters_list <- vector(mode = "list", length = number_of_params)
  
  for (o in 1:length(partial_states_list)) {
    
    for (p in 1:length(filters_list)) {
      
      filters_list[[p]] <- which(simulated_filter_state_params[,p] >= partial_states_list[[o]][[p]][[1]] & simulated_filter_state_params[,p] <= partial_states_list[[o]][[p]][[2]])
      
    }
    
    partial_states_filters_list[[o]] <- filters_list
    
  }
  
  #find intersection of the filters
  filtered_intersections_list <- vector(mode = "list", length = length(partial_states_list))
  percent_models <- vector(length = length(partial_states_list))
  
  for (q in 1:length(partial_states_filters_list)) {
    
    filtered_intersection <- partial_states_filters_list[[q]][[1]]     
    
    for (r in 2:length(partial_states_filters_list[[q]])) {
      
      filtered_intersection <- intersect(filtered_intersection, partial_states_filters_list[[q]][[r]])
      
    }
    
    filtered_intersections_list[[q]] <- filtered_intersection
    percent_models[q] <- 100*length(filtered_intersection)/numModels_filter_state
    
  }
  
  #filter the filter_state and calculate score
  partial_state_scores <- vector(mode = "list", length = length(partial_states_list))
  partial_state_scores_sums <- vector(length = length(partial_states_list))
  
  for (s in 1:length(partial_state_scores)) {
    
    filtered_expr_log <- log2(t(assay(simulated_filter_state$simulated_rset)))[filtered_intersections_list[[s]],]
    
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
      
      #ks score for distribution of each gene
      for (w in 1:length(cf_vec)) {
        
        ks_test <- dgof::ks.test(x = reference_rset_object$prcomp_ref$x[,w], y = filtered_expr_PCs[,w])
        cf_vec[w] <- as.numeric(ks_test$statistic)
        
      }
      
    } else {
      
      for (w in 1:length(cf_vec)) {
        
        ks_test <- dgof::ks.test(x = reference_rset_object$expression_ref[,w], y = filtered_expr_log_z[,w])
        cf_vec[w] <- as.numeric(ks_test$statistic)
        
      }
    }
    
    partial_state_scores[[s]] <- cf_vec  
    partial_state_scores_sums[s] <- sum(cf_vec)
    
    
  }
  
  if (!missing(current_state_costs_object)) {
    
    current_state_costs_object <- current_state_costs_object
  } else {
    
  #calculate current state score
  current_state_costs_object <- calc_score(state = current_state, reference_rset_object = reference_rset_object, scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, save_rset = TRUE)
  }
  
  #calculate derivatives
  #subtract partial state score from current state score
  
  partial_derivs_vec <- vector(length = length(partial_states_list))
  
  for (a in 1:length(partial_derivs_vec)) {
    
    partial_derivs_vec[a] <- (partial_state_scores_sums[a] - current_state_costs_object$cf_vec_sum)/update_vector_rep[a]
    
  }
  
  if (full_save == TRUE) {
    
    partial_derivs_object <- list(update_vector = update_vector, current_state_costs_object = current_state_costs_object, partial_states_list = partial_states_list, simulated_filter_state = simulated_filter_state, 
                                  partial_states_filters_list = partial_states_filters_list, filtered_intersections_list = filtered_intersections_list, 
                                  partial_state_scores = partial_state_scores, partial_state_scores_sums = partial_state_scores_sums, partial_derivs_vec = partial_derivs_vec, percent_models = percent_models)
    
  } else {
    

    partial_derivs_object <- list(update_vector = update_vector, partial_state_scores = partial_state_scores, partial_state_scores_sums = partial_state_scores_sums, partial_derivs_vec = partial_derivs_vec, percent_models = percent_models)
        
    
  }

  partial_derivs_object
  
}

#####################partial_derivs_loop##################################################################################################################################################################
partial_derivs_loop <- function(number_of_iterations = 250, current_state, reference_rset_object, global_bounds, step_size = 0.01, numModels = 10000, integrateStepSize = 0.02, numModels_filter_state = 20000, integrateStepSize_filter_state = 0.02, scoring_method = scoring_method,
                                use_global_bounds = FALSE, full_save = FALSE, resimulate_current_state = FALSE) {
  
  
  if (resimulate_current_state == FALSE) {
    
    current_state_costs_object <- calc_score(state = current_state, reference_rset_object = reference_rset_object, scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, save_rset = TRUE)
    
    partial_derivs_loop_list <- vector(mode = "list", length = number_of_iterations)
    
    for (i in 1:number_of_iterations) {
      
      partial_derivs_loop_list[[i]] <- calc_partial_derivs(current_state = current_state, reference_rset_object = reference_rset_object, global_bounds = global_bounds, step_size = step_size, numModels = numModels, integrateStepSize = integrateStepSize, 
                                                             numModels_filter_state = numModels_filter_state, integrateStepSize_filter_state = integrateStepSize_filter_state, scoring_method = scoring_method, use_global_bounds = use_global_bounds, full_save = full_save, 
                                                             current_state_costs_object = current_state_costs_object)
      
    }
    
    
    
  } else {
  
  partial_derivs_loop_list <- vector(mode = "list", length = number_of_iterations)
  
  for (i in 1:number_of_iterations) {
    
    partial_derivs_loop_list[[i]] <- calc_partial_derivs(current_state = current_state, reference_rset_object = reference_rset_object, global_bounds = global_bounds, step_size = step_size, numModels = numModels, integrateStepSize = integrateStepSize, 
                                                           numModels_filter_state = numModels_filter_state, integrateStepSize_filter_state = integrateStepSize_filter_state, scoring_method = scoring_method, use_global_bounds = use_global_bounds, full_save = full_save)
    
  }
  }
  partial_derivs_loop_object <- list(partial_derivs_loop_list = partial_derivs_loop_list, arguments = list(current_state = current_state, reference_rset_object = reference_rset_object, global_bounds = global_bounds,
                                                                                                           step_size = step_size, numModels = numModels, integrateStepSize = integrateStepSize,
                                                                                                           numModels_filter_state = numModels_filter_state, integrateStepSize_filter_state = integrateStepSize_filter_state, scoring_method = scoring_method,
                                                                                                           use_global_bounds = FALSE, resimulate_current_state = FALSE))
  
}
#####################tests##################################################################################################################################################################
partial_derivs_loop_01_20000_Default_TRUE <- partial_derivs_loop(number_of_iterations = 250, current_state = ref_current_state0, reference_rset_object = test, global_bounds = global_bounds_list, step_size = 0.01, numModels = 10000,
                                            integrateStepSize = 0.02, numModels_filter_state = 20000, integrateStepSize_filter_state = 0.02, scoring_method = "Default",
                                            use_global_bounds = FALSE, full_save = FALSE, resimulate_current_state = TRUE)

partial_derivs_loop_01_10000_Default_TRUE <- partial_derivs_loop(number_of_iterations = 250, current_state = ref_current_state0, reference_rset_object = test, global_bounds = global_bounds_list, step_size = 0.01, numModels = 10000,
                                                                 integrateStepSize = 0.02, numModels_filter_state = 10000, integrateStepSize_filter_state = 0.02, scoring_method = "Default",
                                                                 use_global_bounds = FALSE, full_save = FALSE, resimulate_current_state = TRUE)
  
partial_derivs_loop_03_20000_Default_TRUE <- partial_derivs_loop(number_of_iterations = 250, current_state = ref_current_state0, reference_rset_object = test, global_bounds = global_bounds_list, step_size = 0.03, numModels = 10000,
                                                                 integrateStepSize = 0.02, numModels_filter_state = 20000, integrateStepSize_filter_state = 0.02, scoring_method = "Default",
                                                                 use_global_bounds = FALSE, full_save = FALSE, resimulate_current_state = TRUE)

##################################################################################################################################################################



#tests
a <- calc_partial_derivs(current_state = ref_current_state0, reference_rset_object = test, global_bounds = global_bounds_list, step_size = 0.01, numModels = 10000, integrateStepSize = 0.02, 
                         numModels_filter_state = 10000, integrateStepSize_filter_state = 0.02, scoring_method = "Default")
a$partial_derivs_vec

a2 <- calc_partial_derivs(current_state = ref_current_state0, reference_rset_object = test, global_bounds = global_bounds_list, step_size = 0.01, numModels = 10000, integrateStepSize = 0.02, 
                         numModels_filter_state = 10000, integrateStepSize_filter_state = 0.02, scoring_method = "Default", use_global_bounds = FALSE, full_save = FALSE)

partial_derivs_loop_01_20000_Default_TRUE$partial_derivs_loop_list[[1]]$partial_derivs_vec

partialvec_list <- vector(mode = "list",length = 40)
partialvec <- vector(length = 250)
partial_means <- vector(length = 40)

for (i in 1:40) {
  
  for (j in 1:250) {
    
    partialvec[j] <- partial_derivs_loop_03_20000_Default_TRUE$partial_derivs_loop_list[[j]]$partial_derivs_vec[i]
    
  }


  partialvec_list[[i]] <- partialvec
  partial_means[i] <- mean(partialvec)
}

which(partial_means == -0.2140001183)

sort(partial_means)
order(partial_means)

partialvec1 <- vector(length = 250)
for (k in 1:250) {
  
partialvec  
  
  
  
}

length(partialvec_list)

mean(partialv1)

plot(density(partialvec_list[[26]]))
sd(partialvec_list[[11]])

a2$percent_models

a$update_vector
unlist(a$partial_states_list[[2]])
unlist(a$filter_state)
head(a$simulated_filter_state_params)

i <- 1
vec <- vector(length = length(a$filter_state))
for (i in 1:length(vec)) {
  
  vec[i] <-  a$filter_state[[i]][[2]] - a$filter_state[[i]][[1]]
  
  if (i == 7 | i == 8) {
    
    vec[i] <-  round(a$filter_state[[i]][[2]]) - round(a$filter_state[[i]][[1]]) 
  }
}
vec

i <- 1
vec2 <- vector(length = length(a$filter_state))
for (i in 1:length(vec2)) {
  vec2[i] <-  a$partial_states_list[[4]][[i]][[2]] - a$partial_states_list[[4]][[i]][[1]]
  
  if (i == 7 | i == 8) {
    
    vec2[i] <-  round(a$partial_states_list[[4]][[i]][[2]]) - round(a$partial_states_list[[4]][[i]][[1]])
    
  }
}
vec2
vec2/vec
prod(vec2/vec)*10000
length(a$filtered_intersections_list[[4]])


#Current state: Pa: 20:40, Pb:20,40 (simulate 10k)
#add to low bound, 21 If 21 is in the range of current state, then filter.
#subract from lower bound: 19, 19:40, depends on the relative percent, then combine
#(if lower bound becomes less than global lower bound, then set equal to global lower bound)
#Filter state: Pa


