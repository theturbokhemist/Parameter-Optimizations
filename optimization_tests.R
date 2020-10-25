
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
  if (scoring_method1 == "PCs") {
    
    simulated_expresion <- simulated_rset_expr_PCs
    
    #ks score for distribution of each gene
    for (w in 1:length(cf_vec)) {
      
      ks_test <- dgof::ks.test(x = reference_rset_object$prcomp_ref$x[,w], y = simulated_rset_expr_PCs[,w])
      cf_vec[w] <- as.numeric(ks_test$statistic)
      
    }
    
  } else {
    
    simulated_expresion <- simulated_rset_expr_log_z
    
    for (w in 1:length(cf_vec)) {
      
      ks_test <- dgof::ks.test(x = reference_rset_object$expression_ref[,w], y = simulated_rset_expr_log_z[,w])
      cf_vec[w] <- as.numeric(ks_test$statistic)
      
    }
  }
  
  costs_object <- list(cf_vec = cf_vec, cf_vec_sum = sum(cf_vec))
  
  if (save_rset == TRUE) {
    
    costs_object <- list(cf_vec = cf_vec, cf_vec_sum = sum(cf_vec), simulated_rset = simulated_rset,  simulated_expression = simulated_expression,
                         arguments = list(scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize))
    
  }
  costs_object
}

# calc_score(state = ref_state0, reference_rset_object = test, scoring_method = "PCs", numModels = 10000, integrateStepSize = 0.02)
##################################################################################################################################################################

ref_state0 <- list(list(1,100), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))
ref_current_state0 <- list(list(20,50), list(20,50), list(0.2,0.5), list(0.2, 0.5), list(10,40), list(10,40), list(2,4), list(2,4), list(20,50), list(20,50))

topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))
scoring_method1 <- "Default"
test <- generate_ref_object(topology = topology_TS, reference_state = ref_state0, numModels = 2000, integrateStepSize = 0.02)
global_bounds_list <- list(G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0, TH_upper = 55, N_lower = 1, N_upper = 5, FC_lower = 1, FC_upper = 100)

#####################calc_partial_derivs##################################################################################################################################################################
#1  current state: P_A 10-20, P_B: 10-20, calculate E0

#2: new state: substract step size from every lower bound and add step size to every upper bound

#3 simulate new state: P_A 9-21, P_B: 9-21

#4 create list of all partial states
# a.) P_A 9-20, P_B: 10-20
# b.) P_A 11-20, P_B: 10-20

# c.) P_A 10-19, P_B: 10-20
# d.) P_A 10-21, P_B: 10-20

# e.) P_A 10-20, P_B: 9-20
# f.) P_A 10-20, P_B: 11-20

# e.) P_A 10-20, P_B: 10-19
# f.) P_A 10-20, P_B: 10-21

#5 filter, calculate scores, subtract scores from E0, create vector of partial derivatives for each partial state. If partial state is identical to current state, partial derivative is 0.

#save hill coefficients as fractions, then round them for RACIPE
a <- calc_partial_derivs(current_state = ref_current_state0, reference_rset_object = test, global_bounds = global_bounds_list, step_size = 0.01)

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

calc_partial_derivs <- function(current_state, reference_rset_object, global_bounds, step_size = 0.01, numModels = 10000, integrateStepSize = 0.02, numModels_filter_state = 20000, integrateStepSize_filter_state = 0.02) {
  
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
  
  #vector of step magnitudes for each parameter. The step magnitude is simply the step_size *(Global upper bound - global lowerbound) 
  update_vector <- vector(length = number_of_params)
  
  for (i in 1:number_of_params) {
    
    if (i %in% production_params) {
    
      update_vector[i] <- (global_bounds$G_upper - global_bounds$G_lower)*step_size
    }
    
    if (i %in% degradation_params) {
      
      update_vector[i] <- (global_bounds$K_upper - global_bounds$K_lower)*step_size
    }
    
    if (i %in% threshold_params) {
      
      update_vector[i] <- (global_bounds$TH_upper - global_bounds$TH_lower)*step_size
    }
    
    if (i %in% hill_params) {
      
      update_vector[i] <- (global_bounds$N_upper - global_bounds$N_lower)*step_size
    }
    
    if (i %in% FC_params) {
      
      update_vector[i] <- (global_bounds$FC_upper - global_bounds$FC_lower)*step_size
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
  
  simulated_filter_state <- sracipeSimulate(reference_rset_object$topology, integrate = FALSE,
                                    numModels = numModels_filter_state, genParams = TRUE,
                                    integrateStepSize = integrateStepSize_filter_state)
  
  #assign vaues from uniform samples distributions to each model of simulated rset.  
  for (m in 1:length(filter_state)) {
    
    sracipeParams(simulated_filter_state)[,m] <- runif(numModels_filter_state, min = filter_state[[m]][[1]], max = filter_state[[m]][[2]])  
    
  }
  
  #sample and round all the hill coefficients
  for (n in 1:length(hill_params)) {
    #you need to subtract 0.5 from lower bound and add 0.5 to upper bound in order to get a uniform distribution
    sracipeParams(simulated_filter_state)[,hill_params[n]] <- round(runif(numModels_filter_state, min= round(filter_state[[hill_params[n]]][[1]]) - 0.5, max= round(filter_state[[hill_params[n]]][[2]]) + 0.5))  
    
  }
  
  
  #simulate filter state
  simulated_filter_state <- sracipeSimulate(simulated_filter_state, integrate = TRUE, genParams = FALSE)
  
  
  simulated_filter_state_params <- sracipeParams(simulated_filter_state)
  
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
  
  for (q in 1:length(partial_states_filters_list)) {
    
    filtered_intersection <- partial_states_filters_list[[q]][[1]]     
    
    for (r in 2:length(partial_states_filters_list[[q]])) {
      
      filtered_intersection <- intersect(filtered_intersection, partial_states_filters_list[[q]][[r]])
      
    }
    
    filtered_intersections_list[[q]] <- filtered_intersection
    
  }
  
  #filter the filter_state and calculate score
  partial_state_scores <- vector(mode = "list", length = length(partial_states_list))
  partial_state_scores_sums <- vector(length = length(partial_states_list))
  
  for (s in 1:length(partial_state_scores)) {
    
    filtered_expr_log <- log2(t(assay(simulated_filter_state)))[filtered_intersections_list[[s]],]
    
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
  
  #calculate derivatives
  #subtract partial state score from current state score
 
  #calculate current state score
  current_state_costs_object <- calc_score(state = current_state, reference_rset_object = reference_rset_object, scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, save_rset = TRUE)
  
  partial_derivs_object <- list(update_vector = update_vector, partial_states_list = partial_states_list, filter_state = filter_state, simulated_filter_state = simulated_filter_state, 
                                simulated_filter_state_params = simulated_filter_state_params, partial_states_filters_list = partial_states_filters_list, filtered_intersections_list = filtered_intersections_list, 
                                partial_state_scores = partial_state_scores, partial_state_scores_sums = partial_state_scores_sums)
  
  
  # 
  # filter <- which(sracipeParams(ref_object_40k$reference_rset)[,1] <= 10)
  # 
  # filtered_expr_log <- log2(t(assay(ref_object_40k$reference_rset)))[filter,]
  # 
  # #z normalize
  # filtered_expr_log_z <- sweep(filtered_expr_log, 
  #                              2, ref_object_40k$means_ref, FUN = "-")
  # 
  # filtered_expr_log_z <- sweep(filtered_expr_log_z, 
  #                              2, ref_object_40k$sds_ref, FUN = "/")
  # 
  # #rotate
  # filtered_expr_PCs <- filtered_expr_log_z %*% ref_object_40k$prcomp_ref$rotation
  # 
  # #cost_function vector
  # cf_vec1 <- vector(length = ncol(filtered_expr_log_z))
  
  
  
}



unlist(ref_state0)

#1  simulated P_A 1-21, P_B: 1-21

#2 filter P_A 1-20, P_B: 1-20, calculate score E_0

#3 P_A 2-20, P_B: 1-20, calculate score E_1

#4 filter P_A 1-21, P_B: 1-20, calculate score E_2

#5 filter P_A 1-20, P_B: 2-20, calculate score E_3

#6 filter P_A 1-20, P_B: 1-21, calculate score E_3


ref_object_40k <- generate_ref_object(topology = topology_TS, reference_state = ref_state0, numModels = 40000, integrateStepSize = 0.02)

filter <- which(sracipeParams(ref_object_40k$reference_rset)[,1] <= 10)

filtered_expr_log <- log2(t(assay(ref_object_40k$reference_rset)))[filter,]

#z normalize
filtered_expr_log_z <- sweep(filtered_expr_log, 
                              2, ref_object_40k$means_ref, FUN = "-")

filtered_expr_log_z <- sweep(filtered_expr_log_z, 
                              2, ref_object_40k$sds_ref, FUN = "/")

#rotate
filtered_expr_PCs <- filtered_expr_log_z %*% ref_object_40k$prcomp_ref$rotation

#cost_function vector
cf_vec1 <- vector(length = ncol(filtered_expr_log_z))

#scoring
if (scoring_method1 == "PCs") {
  
  #ks score for distribution of each gene
  for (l in 1:length(cf_vec1)) {
    
    ks_test <- dgof::ks.test(x = ref_object_40k$prcomp_ref$x[,l], y = filtered_expr_PCs[,l])
    cf_vec1[l] <- as.numeric(ks_test$statistic)
    
  }
  
} else {
  
  for (l in 1:length(cf_vec1)) {
    
    ks_test <- dgof::ks.test(x = ref_object_40k$expression_ref[,l], y = filtered_expr_log_z[,l])
    cf_vec1[l] <- as.numeric(ks_test$statistic)
    
  }
  
  
}
test10 <- sum(cf_vec1)
  

filter <- which(sracipeParams(ref_object_40k$reference_rset)[,1] <= 11)

filtered_expr_log <- log2(t(assay(ref_object_40k$reference_rset)))[filter,]

#z normalize
filtered_expr_log_z <- sweep(filtered_expr_log, 
                             2, ref_object_40k$means_ref, FUN = "-")

filtered_expr_log_z <- sweep(filtered_expr_log_z, 
                             2, ref_object_40k$sds_ref, FUN = "/")

#rotate
filtered_expr_PCs <- filtered_expr_log_z %*% ref_object_40k$prcomp_ref$rotation

#cost_function vector
cf_vec1 <- vector(length = ncol(filtered_expr_log_z))

#scoring
if (scoring_method1 == "PCs") {
  
  #ks score for distribution of each gene
  for (l in 1:length(cf_vec1)) {
    
    ks_test <- dgof::ks.test(x = ref_object_40k$prcomp_ref$x[,l], y = filtered_expr_PCs[,l])
    cf_vec1[l] <- as.numeric(ks_test$statistic)
    
  }
  
} else {
  
  for (l in 1:length(cf_vec1)) {
    
    ks_test <- dgof::ks.test(x = ref_object_40k$expression_ref[,l], y = filtered_expr_log_z[,l])
    cf_vec1[l] <- as.numeric(ks_test$statistic)
    
  }
  
  
}
test11 <- sum(cf_vec1)

test99

test98

test98-test99

  
assay(ref_object_40k$reference_rset)


head(ref_object_40k$expression_ref)
sracipeParams(ref_object_40k$reference_rset)

