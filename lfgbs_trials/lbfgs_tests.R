#####################Globals##################################################################################################################################################################
current_state_og <- list(list(20,50), list(20,50), list(0.2,0.5), list(0.2, 0.5), list(10,40), list(10,40), list(2,4), list(2,4), list(20,50), list(20,50))
current_state_og <- split(unlist(current_state_og), ceiling(seq_along(unlist(current_state_og))/2))

current_state_og2 <- list(list(20,50), list(20,50), list(0.2,0.5), list(0.2, 0.5), list(10,40), list(10,40), list(1.95,4), list(2,4.05), list(20,50), list(20,50))

#######################################################################################################################################################################################
current_state <- list(list(20,50), list(20,50), list(0.2,0.5), list(0.2, 0.5), list(10,40), list(10,40), list(20,50), list(20,50))
state_vector <- unlist(current_state)
skip <- c("N")
#######################################################################################################################################################################################
current_state <- list(list(20,50), list(20,50), list(10,40), list(10,40), list(20,50), list(20,50))
state_vector <- unlist(current_state)
skip <- c("K","N")
#######################################################################################################################################################################################
current_state <- list(list(20,50), list(20,50))
state_vector <- unlist(current_state)
skip <- c("G","K","TH", "N")


#######################################################################################################################################################################################
reference_state <- list(list(1,100), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))
topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))
numModels <- 10000
integrateStepSize <- 0.02
reference <- generate_ref_object(topology = topology_TS, reference_state = reference_state, numModels = numModels, integrateStepSize = integrateStepSize)

global_bounds <- list(G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0, TH_upper = 55, N_lower = 1, N_upper = 5, FC_lower = 1, FC_upper = 100)
g_bounds <- list(G = c(1, 100), K = c(0.1, 1), TH = c(0, 55), N = c(1, 5), FC = c(1, 100))
scoring_method <- "Default"
step_size <- 0.02
gv_numModels_FC <- 5

iteration_pdf <- 0
iteration_cf <- 0

partial_derivs_input <- list()
cost_function_input <- list()

#######################################################################################################################################################################################

test1_lbfgs <- lbfgs::lbfgs(cost_function_wrapper, calculate_partial_derivatives_wrapper, state_vector, min_step = 0.1)
test1_lbfgs$par

test1_lbfgs$par

out.optim <- optim(state_vector, cost_function_wrapper, calculate_partial_derivatives_wrapper, method="L-BFGS-B")
out.optim$par


lbfgs::lbfgs

partial_derivs_input

# install.packages("stochQN")
#######################################################################################################################################################################################
test3 <- stochQN::adaQN(x0 = state_vector, grad_fun = calculate_partial_derivatives_wrapper, obj_fun = cost_function_wrapper)
test3$prev_iter


