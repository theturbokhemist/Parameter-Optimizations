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


#######################################################################################################################################
current_state_below <- list(list(10,40), list(10,40), list(0.15,0.45), list(0.15, 0.45), list(5,20), list(5,20), list(2,3), list(2,3), list(10,40), list(10,40))

ref_state2 <- list(list(50,100), list(50,100), list(0.55,1), list(0.55, 1), list(30,55), list(30,55), list(4,5), list(4,5), list(50,100), list(50,100))

topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))
reference_object2 <- generate_ref_object(topology = topology_TS, reference_state = ref_state2, numModels = 10000, integrateStepSize = 0.02)
global_bounds_list <- list(G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0, TH_upper = 55, N_lower = 1, N_upper = 5, FC_lower = 1, FC_upper = 100)

# partial_derivs_below <- calc_partial_derivatives(current_state = current_state_below, reference_rset_object = reference_object2, global_bounds = global_bounds_list, step_size = 0.001,
#                                                  numModels_filter_state = 10000, integrateStepSize_filter_state = 0.02, scoring_method = "Default", save_all = T)


#################

partial_derivs_calc3 <- function(number_of_iterations, current_state, reference_rset_object, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 500, scoring_method, plots_list, global_bounds, step_size, repeat_current = FALSE) {
  
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
  
  i <- 1
  for (i in 1:number_of_iterations) {
    
    if (repeat_current) {
      
      current_state_scores_object <- calc_score(state = current_state, reference_rset_object = reference_rset_object, scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, save_rset = TRUE)
      current_state_expr_rep <- do.call(rbind, replicate(rep_cstate_number, current_state_scores_object$simulated_expression, simplify = FALSE))
      
    }
    
    
    partial_derivs_vec_L <- c()
    partial_derivs_vec_U <- c()
    
    for (j in 1:length(partial_states_list)) {
      
      partial_state_scores_object <- calc_score(state = partial_states_list[[j]], reference_rset_object = reference_rset_object, scoring_method = scoring_method, numModels = numModels_partial, integrateStepSize = integrateStepSize, save_rset = TRUE, update_vector)
      
      combined_expr <- rbind(current_state_expr_rep, partial_state_scores_object$simulated_expression)
      
      #####
        
        cf_vec <- vector(length = ncol(current_state_scores_object$simulated_expression))
        
        for (w in 1:length(cf_vec)) {
          
          ks_test <- dgof::ks.test(x = reference_rset_object$expression_ref[,w], y = combined_expr[,w])
          cf_vec[w] <- as.numeric(ks_test$statistic)
          cf_sum <-sum(cf_vec)
        }
      
      
      if (j%%2 != 0) {
        
        partial_deriv_minus<- -(cf_sum - current_state_scores_object$cf_vec_sum)/update_vector[j]
        
        partial_derivs_vec_L <- c(partial_derivs_vec_L, partial_deriv_minus)
        
        
      } else {
        
        partial_deriv_plus <- (cf_sum - current_state_scores_object$cf_vec_sum)/update_vector[j]
        
        partial_derivs_vec_U <- c(partial_derivs_vec_U, partial_deriv_plus)
        
      }
      
    }
    
    partial_derivs_list_L[[i]] <- partial_derivs_vec_L
    partial_derivs_list_U[[i]] <- partial_derivs_vec_U
    
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
  
  ##################
  param <- rep(1:10,each = 2)
  
  k <- 1
  
  for (k in 1:length(plots_list)) {
    
    if (k %%2 != 0) {
      
      df <- data.frame(derivs = partial_derivs_list_minus[[param[k]]])
      
      plot <- plots_list[[k]]
      
      plot <- plot + geom_density(df, mapping = aes(x = derivs), color = "green")
      
      plots_list[[k]] <- plot
      
    } else {
      
      df <- data.frame(derivs = partial_derivs_list_plus[[param[k]]])
      
      plot <- plots_list[[k]]
      
      plot <- plot + geom_density(df, mapping = aes(x = derivs), color = "green")
      
      plots_list[[k]] <- plot
      
    }
    
  }
  
  partial_derivs_object <- list(plots_list = plots_list, partial_derivs_L_byiteration = partial_derivs_list_L, partial_derivs_U_byiteration = partial_derivs_list_U, partial_derivs_minus_byparam = partial_derivs_list_minus,
                                partial_derivs_plus_byparam = partial_derivs_list_plus, current_state_scores_object = current_state_scores_object, current_state_expr_rep = current_state_expr_rep, numModels_partial = numModels_partial,
                                rep_cstate_number = rep_cstate_number, partial_states_list = partial_states_list, update_vector = update_vector, combined_expr = combined_expr,
                                arguments = list(current_state = current_state, reference_rset_object = reference_rset_object, numModels = numModels, integrateStepSize = integrateStepSize, numModels_partial = numModels_partial,
                                                 scoring_method = scoring_method, global_bounds = global_bounds_list, step_size = step_size, repeat_current = repeat_current))
  
}





########0.005
partial_derivs_below_ss0005_nMp1000_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 1000,
                                                              scoring_method = "Default", plots_list = stats_test3_0005$plots_list, global_bounds = global_bounds_list, step_size = 0.005, repeat_current = FALSE)

partial_derivs_below_ss0005_nMp500_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 500,
                                                               scoring_method = "Default", plots_list = stats_test3_0005$plots_list, global_bounds = global_bounds_list, step_size = 0.005, repeat_current = FALSE)

partial_derivs_below_ss0005_nMp200_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 200,
                                                               scoring_method = "Default", plots_list = stats_test3_0005$plots_list, global_bounds = global_bounds_list, step_size = 0.005, repeat_current = FALSE)

partial_derivs_below_ss0005_nMp100_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 100,
                                                               scoring_method = "Default", plots_list = stats_test3_0005$plots_list, global_bounds = global_bounds_list, step_size = 0.005, repeat_current = FALSE)


partial_derivs_below_ss0005_nMp50_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 50,
                                                              scoring_method = "Default", plots_list = stats_test3_0005$plots_list, global_bounds = global_bounds_list, step_size = 0.005, repeat_current = FALSE)

list_of_partialderiv_objects_below_ss0005 <- list(partial_derivs_below_ss0005_nMp500_rcF, partial_derivs_below_ss0005_nMp200_rcF, partial_derivs_below_ss0005_nMp1000_rcF, partial_derivs_below_ss0005_nMp100_rcF, partial_derivs_below_ss0005_nMp50_rcF)

test <- extract_sds(list_of_partialderiv_objects = list_of_partialderiv_objects_below_ss0005, reference_rset_object = reference_object2)

do.call("grid.arrange", test$plots_list)

########002

partial_derivs_below_ss002_nMp200_rcT <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 200,
                                                              scoring_method = "Default", plots_list = stats_test3_002$plots_list, global_bounds = global_bounds_list, step_size = 0.02, repeat_current = TRUE)

partial_derivs_below_ss002_nMp200_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 200,
                                                              scoring_method = "Default", plots_list = stats_test3_002$plots_list, global_bounds = global_bounds_list, step_size = 0.02, repeat_current = FALSE)

do.call("grid.arrange", partial_derivs_below_ss002_nMp1000_rcF$plots_list)




partial_derivs_below_ss002_nMp400_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 400,
                                                              scoring_method = "Default", plots_list = stats_test3_002$plots_list, global_bounds = global_bounds_list, step_size = 0.02, repeat_current = FALSE)

partial_derivs_below_ss002_nMp600_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 600,
                                                              scoring_method = "Default", plots_list = stats_test3_002$plots_list, global_bounds = global_bounds_list, step_size = 0.02, repeat_current = FALSE)

partial_derivs_below_ss002_nMp1000_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 1000,
                                                              scoring_method = "Default", plots_list = stats_test3_002$plots_list, global_bounds = global_bounds_list, step_size = 0.02, repeat_current = FALSE)

partial_derivs_below_ss002_nMp2000_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 2000,
                                                               scoring_method = "Default", plots_list = stats_test3_002$plots_list, global_bounds = global_bounds_list, step_size = 0.02, repeat_current = FALSE)

list_of_partialderiv_objects_below_ss002 <- list(partial_derivs_below_ss002_nMp200_rcF, partial_derivs_below_ss002_nMp400_rcF, partial_derivs_below_ss002_nMp600_rcF, partial_derivs_below_ss002_nMp2000_rcF)

test <- extract_sds(list_of_partialderiv_objects = list_of_partialderiv_objects_below_ss002, reference_rset_object = reference_object2)

do.call("grid.arrange", test$plots_list)


########005


partial_derivs_below_ss005_nMp500_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 500,
                                                               scoring_method = "Default", plots_list = stats_test3$plots_list, global_bounds = global_bounds_list, step_size = 0.05, repeat_current = FALSE)

partial_derivs_below_ss005_nMp1000_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 1000,
                                                              scoring_method = "Default", plots_list = stats_test3$plots_list, global_bounds = global_bounds_list, step_size = 0.05, repeat_current = FALSE)

partial_derivs_below_ss005_nMp1500_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 1500,
                                                              scoring_method = "Default", plots_list = stats_test3$plots_list, global_bounds = global_bounds_list, step_size = 0.05, repeat_current = FALSE)

partial_derivs_below_ss005_nMp2500_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 2500,
                                                              scoring_method = "Default", plots_list = stats_test3$plots_list, global_bounds = global_bounds_list, step_size = 0.05, repeat_current = FALSE)

partial_derivs_below_ss005_nMp5000_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 5000,
                                                              scoring_method = "Default", plots_list = stats_test3$plots_list, global_bounds = global_bounds_list, step_size = 0.05, repeat_current = FALSE)

list_of_partialderiv_objects_below_ss005 <- list(partial_derivs_below_ss005_nMp500_rcF, partial_derivs_below_ss005_nMp1000_rcF, partial_derivs_below_ss005_nMp1500_rcF, partial_derivs_below_ss005_nMp2500_rcF)

########001


partial_derivs_below_ss001_nMp100_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 100,
                                                              scoring_method = "Default", plots_list = stats_test3_001$plots_list, global_bounds = global_bounds_list, step_size = 0.01, repeat_current = FALSE)

partial_derivs_below_ss001_nMp200_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 200,
                                                               scoring_method = "Default", plots_list = stats_test3_001$plots_list, global_bounds = global_bounds_list, step_size = 0.01, repeat_current = FALSE)

partial_derivs_below_ss001_nMp300_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 300,
                                                               scoring_method = "Default", plots_list = stats_test3_001$plots_list, global_bounds = global_bounds_list, step_size = 0.01, repeat_current = FALSE)

partial_derivs_below_ss001_nMp500_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 500,
                                                               scoring_method = "Default", plots_list = stats_test3_001$plots_list, global_bounds = global_bounds_list, step_size = 0.01, repeat_current = FALSE)

partial_derivs_below_ss001_nMp1000_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 1000,
                                                               scoring_method = "Default", plots_list = stats_test3_001$plots_list, global_bounds = global_bounds_list, step_size = 0.01, repeat_current = FALSE)

list_of_partialderiv_objects_below_ss001 <- list(partial_derivs_below_ss001_nMp100_rcF, partial_derivs_below_ss001_nMp200_rcF, partial_derivs_below_ss001_nMp300_rcF, partial_derivs_below_ss001_nMp500_rcF, partial_derivs_below_ss001_nMp1000_rcF)

########0001


partial_derivs_below_ss0001_nMp10_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 10,
                                                              scoring_method = "Default", plots_list = stats_test3_0001$plots_list, global_bounds = global_bounds_list, step_size = 0.001, repeat_current = FALSE)

partial_derivs_below_ss0001_nMp20_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 20,
                                                              scoring_method = "Default", plots_list = stats_test3_0001$plots_list, global_bounds = global_bounds_list, step_size = 0.001, repeat_current = FALSE)

partial_derivs_below_ss0001_nMp30_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 30,
                                                              scoring_method = "Default", plots_list = stats_test3_0001$plots_list, global_bounds = global_bounds_list, step_size = 0.001, repeat_current = FALSE)

partial_derivs_below_ss0001_nMp50_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 50,
                                                              scoring_method = "Default", plots_list = stats_test3_0001$plots_list, global_bounds = global_bounds_list, step_size = 0.001, repeat_current = FALSE)

partial_derivs_below_ss0001_nMp100_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 100,
                                                               scoring_method = "Default", plots_list = stats_test3_0001$plots_list, global_bounds = global_bounds_list, step_size = 0.001, repeat_current = FALSE)

partial_derivs_below_ss0001_nMp500_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_below, reference_rset_object = reference_object2, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 500,
                                                               scoring_method = "Default", plots_list = stats_test3_0001$plots_list, global_bounds = global_bounds_list, step_size = 0.001, repeat_current = FALSE)

list_of_partialderiv_objects_below_ss0001 <- list(partial_derivs_below_ss0001_nMp10_rcF, partial_derivs_below_ss0001_nMp20_rcF, partial_derivs_below_ss0001_nMp30_rcF, partial_derivs_below_ss0001_nMp50_rcF,
                                                  partial_derivs_below_ss0001_nMp100_rcF, partial_derivs_below_ss0001_nMp500_rcF)

do.call("grid.arrange", partial_derivs_below_ss005_nMp2500_rcF$plots_list)

list_of_partialderiv_objects_below_ss005 <- list(partial_derivs_below_ss005_nMp500_rcF, partial_derivs_below_ss005_nMp1000_rcF, partial_derivs_below_ss005_nMp1500_rcF, partial_derivs_below_ss005_nMp2500_rcF)

test1 <- extract_sds(list_of_partialderiv_objects = list_of_partialderiv_objects_below_ss005, reference_rset_object = reference_object2)

do.call("grid.arrange", test$plots_list)

test$plots_list[[1]] + test1$plots_list[[1]]


do.call("grid.arrange", partial_derivs_below_ss001_nMp500_rcF$plots_list)

# 
# partial_derivs_below_ss001_nMp1000_rcF <- partial_derivs_below_ss0005_nMp1000_rcF
# 
# partial_derivs_below_ss001_nMp500_rcF <- partial_derivs_below_ss0005_nMp500_rcF
# 
# partial_derivs_below_ss001_nMp200_rcF <- partial_derivs_below_ss0005_nMp200_rcF
# 
# partial_derivs_below_ss001_nMp100_rcF <- partial_derivs_below_ss0005_nMp100_rcF

partial_derivs_below_ss001_nMp100_rcF$arguments$

do.call("grid.arrange", partial_derivs_below_ss0005_nMp1000_rcF$plots_list)


extract_sds <- function(list_of_partialderiv_objects, reference_rset_object, fold_change = FALSE, relativeCV = TRUE, divide_by_all_means = TRUE) {
  
  #vector of paramater names
  parameter_names <- colnames(sracipeParams(reference_rset_object$reference_rset))
  
  #vector of parameter_names + bound
  parameter_bounds_names <- paste(rep(parameter_names, each = 2),
                                  rep(c("lower", "upper"), times = length(parameter_names)), sep = "_")
  
  
  step_size <- list_of_partialderiv_objects[[1]]$arguments$step_size
  
  ##
  
  addedModels <- c()
  
  for (i in 1:length(list_of_partialderiv_objects)) {
    
    if (fold_change) {
      
      addedModels <- c(addedModels, list_of_partialderiv_objects[[i]]$arguments$numModels_partial/((list_of_partialderiv_objects[[i]]$arguments$step_size)*(list_of_partialderiv_objects[[i]]$arguments$numModels)))
      
    } else {
    
      addedModels <- c(addedModels, list_of_partialderiv_objects[[i]]$arguments$numModels_partial)
    
    }
    
  }
  
  
  list_of_partialderiv_objects <- list_of_partialderiv_objects[order(addedModels)]
  
  addedModels <- sort(addedModels)
  
  df <- data.frame(addedModels = addedModels, percentCV = 1:length(list_of_partialderiv_objects))
  
  df_means <- data.frame(addedModels = addedModels, means = 1:length(list_of_partialderiv_objects))
  
  plots_list <- vector(mode = "list", length = 20)
  
  means_list <- vector(mode = "list", length = 20)
  
  dfs_list <- vector(mode = "list", length = 20)
  
  dfs_list_means <- vector(mode = "list", length = 20)
  
  ###
  
  param <- rep(1:10,each = 2)
  
  for (j in 1:20) {
    
    means <- vector(length = length(list_of_partialderiv_objects))
    
    if (j %%2 != 0) {
      
      ###
      
      for (k in 1:length(list_of_partialderiv_objects)) {
        
        means[k] <-abs(mean(list_of_partialderiv_objects[[k]]$partial_derivs_minus_byparam[[param[j]]]))
        
      }
      
      df_means[, 2] <- means
      
      
      ###
      
      for (k in 1:length(list_of_partialderiv_objects)) {
        
        if (divide_by_all_means) {
          
          df[k, 2] <- (sd(list_of_partialderiv_objects[[k]]$partial_derivs_minus_byparam[[param[j]]]))/mean(means)
          
        } else {
          
          df[k, 2] <- (sd(list_of_partialderiv_objects[[k]]$partial_derivs_minus_byparam[[param[j]]]))/means[k]
          
        }
        

      }
      
      ###
      
      if (relativeCV) {
        
        df[, 2] <- df[, 2]/max(df[, 2])

      }
      
      plot <- ggplot(data = df, aes(x = addedModels, y = percentCV)) + geom_line() + geom_point() + labs(title=parameter_bounds_names[j])
      
      
    } else {
      
      ###
      
      for (k in 1:length(list_of_partialderiv_objects)) {
        
        means[k] <-abs(mean(list_of_partialderiv_objects[[k]]$partial_derivs_plus_byparam[[param[j]]]))
        
      }
      
      df_means[, 2] <- means
      ###
      
      for (k in 1:length(list_of_partialderiv_objects)) {
        
        if (divide_by_all_means) {
        
        df[k, 2] <- (sd(list_of_partialderiv_objects[[k]]$partial_derivs_plus_byparam[[param[j]]]))/mean(means)
        
        } else {
          
          df[k, 2] <- (sd(list_of_partialderiv_objects[[k]]$partial_derivs_plus_byparam[[param[j]]]))/means[k]
          
        }
        
      }
      
      ###
      
      if (relativeCV) {
        
        df[, 2] <- df[, 2]/max(df[, 2])
        
      }
      
      
      plot <- ggplot(data = df, aes(x = addedModels, y = percentCV)) + geom_line() + geom_point() + labs(title=parameter_bounds_names[j])
      
    }

    dfs_list_means[[j]] <- df_means
    dfs_list[[j]] <- df
    plots_list[[j]] <- plot
    means_list[[j]] <- means
    
  }
  object <- list(plots_list = plots_list, means_list = means_list, dfs_list = dfs_list, step_size = step_size, dfs_list_means = dfs_list_means, addedModels = df[,1])
  
}


list_of_partialderiv_objects_below_ss0005 <- list(partial_derivs_below_ss0005_nMp500_rcF, partial_derivs_below_ss0005_nMp200_rcF, partial_derivs_below_ss0005_nMp1000_rcF, partial_derivs_below_ss0005_nMp100_rcF)

test <- extract_sds(list_of_partialderiv_objects = list_of_partialderiv_objects_below_ss0005, reference_rset_object = reference_object2, fold_change = T, divide_by_all_means = T)

test$step_size

overlap_plots <-function(list_of_stats_objects, reference_rset_object, plot_means = FALSE) {
  
  #vector of paramater names
  parameter_names <- colnames(sracipeParams(reference_rset_object$reference_rset))
  
  #vector of parameter_names + bound
  parameter_bounds_names <- paste(rep(parameter_names, each = 2),
                                  rep(c("lower", "upper"), times = length(parameter_names)), sep = "_")
  
  ###
  
  plots_list <- vector(mode = "list", length = 20)
  
  colors <- c("blue", "red", "green", "orange", "purple")
  
  for (i in 1:length(plots_list)) {
    
    if (plot_means) {
      
      df <- list_of_stats_objects[[1]]$dfs_list_means[[i]]
      
      
      plot <- ggplot(data = df, aes(x = addedModels, y = means)) + geom_line(color = colors[1] ) + geom_point() + labs(title=parameter_bounds_names[i])
      
    } else {
    
    df <- list_of_stats_objects[[1]]$dfs_list[[i]]
    
    plot <- ggplot(data = df, aes(x = addedModels, y = percentCV)) + geom_line(color = colors[1] ) + geom_point() + labs(title=parameter_bounds_names[i])
    
    }
    
    plots_list[[i]] <- plot
  }
  
  for (j in 2:length(list_of_stats_objects)) {
    
    for (k in 1:length(plots_list)) {
      
      if (plot_means) {
        
        df <- list_of_stats_objects[[j]]$dfs_list_means[[k]]
        
        plot <- plot <- plots_list[[k]] + geom_line(data = df, aes(x = addedModels, y = means) ,color = colors[j] ) + geom_point(data = df, aes(x = addedModels, y = means))
        
      } else {
      
      df <- list_of_stats_objects[[j]]$dfs_list[[k]]
      
      # print(colnames(list_of_stats_objects[[1]]$dfs_list[[i]])[1])
      # 
      # print(colnames(list_of_stats_objects[[1]]$dfs_list[[i]])[2])
      
      plot <- plots_list[[k]] + geom_line(data = df, aes(x = addedModels, y = percentCV) ,color = colors[j] ) + geom_point(data = df, aes(x = addedModels, y = percentCV))
      
      }
      plots_list[[k]] <- plot
    }
    
    
    
  }
  

  plots_list
  
}


CVplot__below_ss0001_fcT <- extract_sds(list_of_partialderiv_objects = list_of_partialderiv_objects_below_ss0001, reference_rset_object = reference_object2, fold_change = T)
CVplot__below_ss0005_fcT <- extract_sds(list_of_partialderiv_objects = list_of_partialderiv_objects_below_ss0005, reference_rset_object = reference_object2, fold_change = T)
CVplot__below_ss001_fcT <- extract_sds(list_of_partialderiv_objects = list_of_partialderiv_objects_below_ss001, reference_rset_object = reference_object2, fold_change = T)
CVplot__below_ss002_fcT <- extract_sds(list_of_partialderiv_objects = list_of_partialderiv_objects_below_ss002, reference_rset_object = reference_object2, fold_change = T)
CVplot__below_ss005_fcT <- extract_sds(list_of_partialderiv_objects = list_of_partialderiv_objects_below_ss005, reference_rset_object = reference_object2, fold_change = T)

overlap_CVplot_below <- overlap_plots(list(CVplot__below_ss0001_fcT, CVplot__below_ss0005_fcT, CVplot__below_ss001_fcT, CVplot__below_ss002_fcT, CVplot__below_ss005_fcT), reference_rset_object = reference_object2, plot_means = TRUE)
means_barplot_below <- means_plot(list(CVplot__below_ss0005_fcT, CVplot__below_ss0001_fcT, CVplot__below_ss001_fcT, CVplot__below_ss002_fcT, CVplot__below_ss005_fcT), reference_rset_object = reference_object2)

test <- means_plot(list(CVplot__below_ss0005_fcT, CVplot__below_ss0001_fcT, CVplot__below_ss001_fcT, CVplot__below_ss002_fcT, CVplot__below_ss005_fcT), reference_rset_object = reference_object2, plot_mean_variation = F, plot_vs_stepsize = F)
test <- overlap_plots(list(CVplot__below_ss0005_fcT, CVplot__below_ss0001_fcT, CVplot__below_ss001_fcT, CVplot__below_ss002_fcT, CVplot__below_ss005_fcT), reference_rset_object = reference_object2, plot_means = T)
do.call("grid.arrange", test)


do.call("grid.arrange", means_barplot_below)

do.call("grid.arrange", overlap_CVplot_below)
test$means_list
#step 1: simulate P_A 10:40, 10k times
#step 2: simulate P_A 9.97: 10 (0.1%), 500 times
#step 3: repeat 10k model expression matrix 50 times to make 500k long expression matrix
#step 4: combine 500k matrix with 500 model matrix
#step 5: ks score between 500.5 k matrix and 10 k matrix


#######cs_between########
########0001

partial_derivs_between_ss0001_nMp10_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 10,
                                                              scoring_method = "Default", plots_list = statistics_csbetween_ss0001$plots_list, global_bounds = global_bounds_list, step_size = 0.001, repeat_current = FALSE)

partial_derivs_between_ss0001_nMp20_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 20,
                                                              scoring_method = "Default", plots_list = statistics_csbetween_ss0001$plots_list, global_bounds = global_bounds_list, step_size = 0.001, repeat_current = FALSE)

partial_derivs_between_ss0001_nMp30_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 30,
                                                              scoring_method = "Default", plots_list = statistics_csbetween_ss0001$plots_list, global_bounds = global_bounds_list, step_size = 0.001, repeat_current = FALSE)

partial_derivs_between_ss0001_nMp50_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 50,
                                                              scoring_method = "Default", plots_list = statistics_csbetween_ss0001$plots_list, global_bounds = global_bounds_list, step_size = 0.001, repeat_current = FALSE)

partial_derivs_between_ss0001_nMp100_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 100,
                                                               scoring_method = "Default", plots_list = statistics_csbetween_ss0001$plots_list, global_bounds = global_bounds_list, step_size = 0.001, repeat_current = FALSE)

partial_derivs_between_ss0001_nMp200_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 200,
                                                                 scoring_method = "Default", plots_list = statistics_csbetween_ss0001$plots_list, global_bounds = global_bounds_list, step_size = 0.001, repeat_current = FALSE)

partial_derivs_between_ss0001_nMp500_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 500,
                                                               scoring_method = "Default", plots_list = statistics_csbetween_ss0001$plots_list, global_bounds = global_bounds_list, step_size = 0.001, repeat_current = FALSE)

list_of_partialderiv_objects_between_ss0001 <- list(partial_derivs_between_ss0001_nMp10_rcF, partial_derivs_between_ss0001_nMp20_rcF, partial_derivs_between_ss0001_nMp30_rcF, partial_derivs_between_ss0001_nMp50_rcF,
                                                    partial_derivs_between_ss0001_nMp100_rcF, partial_derivs_between_ss0001_nMp500_rcF)

########0005

partial_derivs_between_ss0005_nMp50_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 50,
                                                                scoring_method = "Default", plots_list = statistics_csbetween_ss0005$plots_list, global_bounds = global_bounds_list, step_size = 0.005, repeat_current = FALSE)

partial_derivs_between_ss0005_nMp100_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 100,
                                                                scoring_method = "Default", plots_list = statistics_csbetween_ss0005$plots_list, global_bounds = global_bounds_list, step_size = 0.005, repeat_current = FALSE)

partial_derivs_between_ss0005_nMp150_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 150,
                                                                scoring_method = "Default", plots_list = statistics_csbetween_ss0005$plots_list, global_bounds = global_bounds_list, step_size = 0.005, repeat_current = FALSE)

partial_derivs_between_ss0005_nMp250_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 250,
                                                                scoring_method = "Default", plots_list = statistics_csbetween_ss0005$plots_list, global_bounds = global_bounds_list, step_size = 0.005, repeat_current = FALSE)

partial_derivs_between_ss0005_nMp500_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 500,
                                                                 scoring_method = "Default", plots_list = statistics_csbetween_ss0005$plots_list, global_bounds = global_bounds_list, step_size = 0.005, repeat_current = FALSE)

partial_derivs_between_ss0005_nMp1000_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 1000,
                                                                 scoring_method = "Default", plots_list = statistics_csbetween_ss0005$plots_list, global_bounds = global_bounds_list, step_size = 0.005, repeat_current = FALSE)

list_of_partialderiv_objects_between_ss0005 <- list(partial_derivs_between_ss0005_nMp50_rcF, partial_derivs_between_ss0005_nMp100_rcF, partial_derivs_between_ss0005_nMp150_rcF, partial_derivs_between_ss0005_nMp250_rcF,
                                                    partial_derivs_between_ss0005_nMp500_rcF, partial_derivs_between_ss0005_nMp1000_rcF)


########001

partial_derivs_between_ss001_nMp100_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 100,
                                                                scoring_method = "Default", plots_list = statistics_csbetween_ss001$plots_list, global_bounds = global_bounds_list, step_size = 0.01, repeat_current = FALSE)



partial_derivs_between_ss001_nMp200_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 200,
                                                                 scoring_method = "Default", plots_list = statistics_csbetween_ss001$plots_list, global_bounds = global_bounds_list, step_size = 0.01, repeat_current = FALSE)

partial_derivs_between_ss001_nMp300_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 300,
                                                                 scoring_method = "Default", plots_list = statistics_csbetween_ss001$plots_list, global_bounds = global_bounds_list, step_size = 0.01, repeat_current = FALSE)

partial_derivs_between_ss001_nMp500_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 500,
                                                                 scoring_method = "Default", plots_list = statistics_csbetween_ss001$plots_list, global_bounds = global_bounds_list, step_size = 0.01, repeat_current = FALSE)

partial_derivs_between_ss001_nMp1000_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 1000,
                                                                 scoring_method = "Default", plots_list = statistics_csbetween_ss001$plots_list, global_bounds = global_bounds_list, step_size = 0.01, repeat_current = FALSE)

partial_derivs_between_ss001_nMp2000_rcF <- partial_derivs_calc3(number_of_iterations = 200, current_state = current_state_between, reference_rset_object = reference_object0, numModels = 10000, integrateStepSize = 0.02, numModels_partial = 2000,
                                                                  scoring_method = "Default", plots_list = statistics_csbetween_ss001$plots_list, global_bounds = global_bounds_list, step_size = 0.01, repeat_current = FALSE)

list_of_partialderiv_objects_between_ss001 <- list(partial_derivs_between_ss001_nMp100_rcF, partial_derivs_between_ss001_nMp200_rcF, partial_derivs_between_ss001_nMp300_rcF, partial_derivs_between_ss001_nMp500_rcF,
                                                   partial_derivs_between_ss001_nMp1000_rcF, partial_derivs_between_ss001_nMp2000_rcF)


CVplot_between_ss0001_fcT <- extract_sds(list_of_partialderiv_objects = list_of_partialderiv_objects_between_ss0001, reference_rset_object = reference_object0, fold_change = T)
CVplot_between_ss0005_fcT <- extract_sds(list_of_partialderiv_objects = list_of_partialderiv_objects_between_ss0005, reference_rset_object = reference_object0, fold_change = T)
CVplot_between_ss001_fcT <- extract_sds(list_of_partialderiv_objects = list_of_partialderiv_objects_between_ss001, reference_rset_object = reference_object0, fold_change = T)

CVplot_between_ss001_fcT$means_list

overlap_CVplot_between <- overlap_plots(list(CVplot_between_ss0001_fcT, CVplot_between_ss0005_fcT, CVplot_between_ss001_fcT), reference_rset_object = reference_object0, plot_means = TRUE)

do.call("grid.arrange", overlap_CVplot_between)

test <- means_plot(list(CVplot_between_ss0001_fcT, CVplot_between_ss0005_fcT, CVplot_between_ss001_fcT), reference_rset_object = reference_object2, plot_mean_variation = T, plot_vs_stepsize = T)
test <- overlap_plots(list(CVplot_between_ss0001_fcT, CVplot_between_ss0005_fcT, CVplot_between_ss001_fcT), reference_rset_object = reference_object2, plot_means = T)
do.call("grid.arrange", test)




#######means_plot##########################################################################################################################################################

means_plot <- function(list_of_stats_objects, reference_rset_object, plot_mean_variation = FALSE, plot_vs_stepsize = TRUE) {
  
  #vector of paramater names
  parameter_names <- colnames(sracipeParams(reference_rset_object$reference_rset))
  
  #vector of parameter_names + bound
  parameter_bounds_names <- paste(rep(parameter_names, each = 2),
                                  rep(c("lower", "upper"), times = length(parameter_names)), sep = "_")
  
  ###
  
  plots_list <- vector(mode = "list", length = 20)
  
  step_sizes <- c()
  
  for (i in 1:length(list_of_stats_objects)) {
    
      
    step_sizes <- c(step_sizes, list_of_stats_objects[[i]]$step_size)

  }
  
  list_of_stats_objects <- list_of_stats_objects[order(step_sizes)]
  
  step_sizes <- sort(step_sizes)
  
  #####
  
  if (plot_mean_variation) {
    
    if (plot_vs_stepsize) {
      
      
      df <- data.frame(step_sizes = as.factor(step_sizes), sds = 1:length(list_of_stats_objects))
      
      
      for (i in 1:length(plots_list)) {
        
        for (j in 1:length(list_of_stats_objects)) {

          
          df[j, 2] <- sd(list_of_stats_objects[[j]]$means_list[[i]])
          
        }

        
        plot <-  ggplot(data = df, aes(x = step_sizes, y = sds, group = 1)) + geom_line() + geom_point() + labs(title=parameter_bounds_names[i])
        
        plots_list[[i]] <- plot
      }
      

    } else {
      
      addedModels <- c()
      
      for (i in 1:length(list_of_stats_objects)) {
        
        addedModels <- c(addedModels, list_of_stats_objects[[i]]$addedModels)
        
      }
      
      addedModels <- sort(unique(addedModels))
      
      addedModels_list <- vector(mode = "list", length = length(addedModels))
      
      
      
      
      for (i in 1:length(addedModels_list)) {
        
        step_size_list <- vector(mode = "list", length = 0)
        
        for (j in 1:length(list_of_stats_objects)) {
          
          means_vec <- c()
          
          
          if (addedModels[i] %in% list_of_stats_objects[[j]]$addedModels) {
            
            
            match <- match(addedModels[i], list_of_stats_objects[[j]]$addedModels)
            
          
            for (k in 1:length(plots_list)) {
              
              means_vec <- c(means_vec, list_of_stats_objects[[j]]$dfs_list_means[[k]][match, 2])
              
            }
          }
          
          if (length(means_vec) > 0) {
            
            step_size_list <- list.append(step_size_list, means_vec)
            
          }
        
          
        }
        
        addedModels_list[[i]] <- step_size_list
        
      }
      
      sd_of_means_list <- vector(mode = "list", length = 0)
      sds_vec <- vector(length = 20)
      addedModels_kept <- c()
      
      for (l in 1:length(addedModels_list)) {
        
        if (length(addedModels_list[[l]]) >= 2) {
          
          addedModels_kept <- c(addedModels_kept, addedModels[l])
          
          for (m in 1:length(plots_list)) {
            
            means_vec <- c()
            
            for (n in 1:length(addedModels_list[[l]])) {
              
              means_vec <- c(means_vec, addedModels_list[[l]][[n]][[m]])
              
            }
            
            sds_vec[m] <- sd(means_vec)
            
          }
          
          sd_of_means_list <- list.append(sd_of_means_list, sds_vec)
          
        }
        
      }
      
      names(sd_of_means_list) <- addedModels_kept
      
      df <- data.frame(addedModels = addedModels_kept, sds_of_means = 1:length(addedModels_kept))
      
      for (i in 1:length(plots_list)) {
        
        for (j in 1:length(sd_of_means_list)) {
          
          df[j,2] <- sd_of_means_list[[j]][[i]]
          
        }
        
        plot <- ggplot(data = df, aes(x = addedModels, y = sds_of_means)) + geom_line() + geom_point() + labs(title=parameter_bounds_names[i])
        
        plots_list[[i]] <- plot
        
      }
      
    }
    
  } else {
  
  
  
  added_models_list <- c()
  
  list_of_stats_objects <- list_of_stats_objects[order(step_sizes)]
  
  step_sizes <- sort(step_sizes)
  
  df <- data.frame(step_sizes = as.factor(step_sizes), means = 1:length(list_of_stats_objects), sds = 1:length(list_of_stats_objects))

  
  for (i in 1:length(plots_list)) {
    
    for (j in 1:length(list_of_stats_objects)) {
      
      df[j, 2] <- mean(list_of_stats_objects[[j]]$means_list[[i]])
      
      df[j, 3] <- sd(list_of_stats_objects[[j]]$means_list[[i]])
      
    }
    
    plot <- ggplot(data=df, aes(x=step_sizes, y=means, fill = step_sizes)) + geom_bar(stat="identity") + geom_errorbar(aes(ymin=means-sds, ymax=means+sds), width=.2) + labs(title=parameter_bounds_names[i])
    
    plots_list[[i]] <- plot
  }
  
  }
  plots_list
}


means_barplot_between <- means_plot(list(CVplot_between_ss0001_fcT, CVplot_between_ss0005_fcT, CVplot_between_ss001_fcT), reference_rset_object = reference_object0)

test <- means_plot(list(CVplot_between_ss0001_fcT, CVplot_between_ss0005_fcT, CVplot_between_ss001_fcT), reference_rset_object = reference_object0, plot_mean_variation = T, plot_vs_stepsize = T)

do.call("grid.arrange", test)

