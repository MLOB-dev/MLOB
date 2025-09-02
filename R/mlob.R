#' Multi-Level Optimal Bayes Function (MLOB)
#'
#'
#' Implements a regularized Bayesian approach that optimizes
#' the estimation of between-group coefficients by minimizing 
#' Mean Squared Error (MSE), balancing both variance and bias.
#' This method provides more reliable estimates 
#' in scenarios with limited data, offering a robust solution for 
#' accurate parameter estimation in multilevel models. The package is designed for researchers 
#' in psychology, education, and related fields who face challenges in 
#' estimating between-group effects in two-level latent variable models, particularly
#' in scenarios with small sample sizes and low intraclass correlation coefficients. 
#'
#' @param formula an object of class "\link{formula}" (or one that can be coerced to that class): a symbolic description of the model to be fitted. Formula specifies the model (e.g., \code{Y ~ X + C...}), where Y is the dependent variable, X is the context variable, which is the focus of most applications of the model  (always included), and C includes all additional covariates.
#' @param data a data frame (or object converted by \link{as.data.frame} to a data frame) containing the variables referenced in the formula. All variables used in the model, including the dependent variable, context variable, covariates, and grouping variable must be present in this data frame.
#' @param group a name of the variable that defines the affiliation of an individual (row) to the specific group.
#' @param balancing.limit a number that represents the threshold of the maximum relative part of the dataset that can be deleted to balance the data. Defaults to \code{0.2}
#' @param conf.level a numeric value representing the confidence level used to calculate confidence intervals for the estimators. Defaults to \code{0.95}, corresponding to a \code{95\%} confidence level.
#' @param jackknife logical variable. If \code{TRUE}, the jackknife re-sampling method will be applied to  calculate the standard error of the between-group and its confidence interval coefficient. Defaults to \code{FALSE}.
#' @param punish.coeff a multiplier that punishes the balancing procedure when deleting the whole group. If punish.coeff is equal to \code{1}, no additional punishment is applied for deleting the group. Higher values intensify the penalty. Defaults to \code{2}.
#' @param ... additional arguments passed to the function.
#'
#' @details
#' This function also verifies whether the data is balanced (i.e., whether each group contains the same number of individuals). If the data is unbalanced, the balancing procedure
#' comes into effect, and identifies the optimal number of items and groups to delete based on the punishment coefficient. If the amount of data deleted is more than defined by threshold
#' (balancing.limit) then results should be interpreted with caution.
#' 
#' The \code{summary()} function produces output similar to:
#' 
#' \preformatted{
#' Summary of Coefficients:
#'                     Estimate Std. Error Lower CI (99%) Upper CI (99%)   Z value   Pr(>|z|) Significance
#' beta_b             0.4279681  0.7544766     -1.5154349       2.371371 0.5672384 0.57055223
#' gamma_Petal.Length 0.4679522  0.2582579     -0.1972762       1.133181 1.8119567 0.06999289            .
#'
#' For comparison, summary of coefficients from unoptimized analysis (ML):
#'                    Estimate   Std. Error Lower CI (99%) Upper CI (99%)      Z value   Pr(>|z|) Significance
#' beta_b             0.6027440 5.424780e+15  -1.397331e+16   1.397331e+16 1.111094e-16 1.00000000
#' gamma_Petal.Length 0.4679522 2.582579e-01  -1.972762e-01   1.133181e+00 1.811957e+00 0.06999289            .
#'
#' Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#' }
#'
#' @return A list containing the results of the regularized Bayesian estimation,
#' which includes the model formula,dependent and context variables, and other relevant details from the analysis.
#' The returned object is of class \code{mlob_result}.
#' 
#' @section Methods:
#' The returned object supports the following S3 methods:
#' \itemize{
#'   \item \code{print(x)} - Display coefficients, standard errors, confidence intervals, Z-values, and p-values
#'   \item \code{summary(x)} - Comprehensive summary with significance stars and comparison to unoptimized ML
#'   \item \code{coef(x)} - Extract coefficients as a data frame
#'   \item \code{se(x)} - Extract standard errors
#'   \item \code{vcov(x)} - Extract variance-covariance matrix (diagonal only)
#'   \item \code{confint(x, parm, level)} - Extract confidence intervals for specified parameters
#'   \item \code{as.data.frame(x)} - Convert results to a data frame format
#'   \item \code{dim(x)} - Return dimensions (number of parameters)
#'   \item \code{length(x)} - Return number of parameters
#'   \item \code{names(x)} - Return parameter names
#'   \item \code{update(x, formula, data, conf.level, balancing.limit, punish.coeff, jackknife)} - Update model with new parameters
#' }
#
#' @author 
#' Valerii Dashuk \email{vadashuk@gmail.com},
#' Binayak Timilsina \email{binayak.timilsina001@gmail.com},
#' Martin Hecht, and
#' Steffen Zitzmann
#' 
#' @references
#' Dashuk, V., Hecht, M., Luedtke, O., Robitzsch, A., & Zitzmann, S. (2024). \doi{10.13140/RG.2.2.18148.39048}
#' 
#' Dashuk, V., Hecht, M., Luedtke, O., Robitzsch, A., & Zitzmann, S. (2024). \doi{10.1007/s41237-025-00264-7}
#' 
#' Luedtke, O., Marsh, H. W., Robitzsch, A., Trautwein, U., Asparouhov, T., & Muthen, B. (2008). \doi{10.1037/a0012869}
#'
#' @examples
#' 
#' # Example 1: usage with the iris dataset
#'
#' result_iris <- mlob(
#' Sepal.Length ~ Sepal.Width + Petal.Length, 
#' data = iris, group = 'Species',
#' conf.level = 0.01,
#' jackknife = FALSE)
#' 
#' # View summary statistics (similar to summary of a linear model);
#' 
#' summary(result_iris)
#' 
#' # Example 2: usage with highly unbalanced mtcars dataset (adjusted balancing.limit)
#' 
#' result_mtcars <- mlob(
#' mpg ~ hp + wt + am + hp:wt + hp:am, 
#' data = mtcars, group = 'cyl', 
#' balancing.limit = 0.35)
#' 
#' # View summary statistics
#' 
#' summary(result_mtcars)
#' 
#' #' # Example 3: Using all available S3 methods on slightly unbalanced ChickWeight dataset
#' 
#' result <- mlob(weight ~ Time, data = ChickWeight, group = 'Diet', jackknife = FALSE)
#' 
#' # Display methods
#' print(result)                    # Display results
#' summary(result)                 # Comprehensive summary
#' coef(result)                    # Extract coefficients
#' se(result)                      # Extract standard errors
#' vcov(result)                    # Extract variance-covariance matrix
#' confint(result)                 # Extract confidence intervals
#' confint(result, "beta_b")       # Extract CI for specific parameter
#' confint(result, level = 0.99)   # Extract CI with different confidence level
#' as.data.frame(result)            # Convert to data frame
#' dim(result)                     # Get dimensions
#' length(result)                  # Get number of parameters
#' names(result)                   # Get parameter names
#' 
#' # Update model with new parameters
#' update(result, conf.level = 0.99)
#' 
#' # List all available methods
#' methods(class = "mlob_result")
#' 
#' @export
mlob <- function(formula, data, group, balancing.limit=0.2, conf.level = 0.95, jackknife = FALSE, punish.coeff = 2, ...) {

  # Save the name of dataframe for displaying the output
  data_name  <- deparse(substitute(data))
  
  # Ensure data is a data frame
  if (!is.data.frame(data)) {
    warning("The 'data' argument is not a data frame. Converting to data frame.")
    data <- as.data.frame(data)
  }

  # Check if all columns in the data are numeric
  #if (!all(sapply(data, is.numeric))) {
   # stop("All columns in the 'data' must be numeric.")
  #}

  # Ensure the 'conf.level' is numeric and between 0 and 1
  if (!is.numeric(conf.level) || length(conf.level) != 1 || conf.level <= 0 || conf.level >= 1) {
    stop("The 'conf.level' argument must be a numeric scalar with value between 0 and 1 (exclusive).")
  }
  
  # Convert confidence level to significance level for internal calculations
  conf.level <- 1 - conf.level

  # Ensure 'jackknife' is a logical (TRUE or FALSE)
  if (!is.logical(jackknife) || length(jackknife) != 1) {
    stop("The 'jackknife' argument must be a scalar Boolean variable, i.e. TRUE or FALSE.")

  }

  #Ensure balancing.limit is a numeric single value between 0 and 1
  if (!is.numeric(balancing.limit) || length(balancing.limit) != 1 || balancing.limit < 0 || balancing.limit > 1) {
    stop("The 'balancing.limit' argument must be a numeric scalar with value between 0 and 1 (inclusive).")
  }

  #Ensure the punish.coeff is numeric single value greater than 0.
  if (!is.numeric(punish.coeff) || length(punish.coeff) != 1 || punish.coeff <= 0) {
    stop("The 'punish.coeff' should be a numeric value greater than 0.")
  }

  
  if (is.factor(data[[group]])) {
    data[[group]] = as.numeric(data[[group]])
  }
  
  # Parse the formula (Y ~ X + C...)
  # to incorporate model.matrix in mlob function,we can modify formula parsing part to automatically handle interaction terms,factors,
  # we parse covariates from the Matrix A not from the data

  all_vars <- all.vars(formula)
  response_var <- all_vars[1]
  predictor_var <- all_vars[2]
  
  # Check that predictor_var is numeric
  if (!is.numeric(data[[predictor_var]])) {
    stop(sprintf("The main predictor variable '%s' must be numeric (not a factor or character).", predictor_var))
  }

  # Create the model matrix using the formula
  A <- stats::model.matrix(formula, data = data)
  A <- as.data.frame(A)

  # Check if intercept exists, rename it to 'y' or create a new column response_var if not present
  if ("(Intercept)" %in% colnames(A)) {
    colnames(A)[colnames(A) == "(Intercept)"] <- response_var
    A[[response_var]] = data[[response_var]]
  } else {
    A[[response_var]] = data[[response_var]] # Creates a new column response_var
    # Move 'response_var' to the first position
    A <- A[, c(response_var, setdiff(colnames(A), response_var))]
  }

  # Get all predictors from the formula (i.e., X + control variables)
  # All columns from the model matrix

  predictors <- colnames(A)

  # Exclude the intercept (if included in the model matrix)
  predictors <- predictors[predictors != response_var]

  # Extract response variable from data(dependent variable)
  response <- data[[response_var]]


  # Ensure the response variable exists in the data
  if (is.null(response)) {
    stop(paste("The response variable", response_var, "is not found in the data."))
  }

  # Now we have `model_matrix` and `response`

  # Extract the group variable from the data using the specified group name
  group_var <- data[[group]]

  # Add the group variable to the model matrix
  A$group <- group_var

  # Check group sizes
  group_counts <- table(A$group)
  # cat("Groups and their sizes:")
  # print(group_counts)

  # Number of groups
  group_num = length(group_counts)


  # Check for balancing
  if (length(unique(A$group)) < 2) {
    stop("Not enough groups to balance the data.")
  }

  # Initialize variables for deletion logic
  tab <- as.numeric(group_counts)  # Convert table to numeric vector
  s <- rep(0, max(tab) - min(tab) + 1)  # Prepare a vector for items to delete
  s2 <- rep(0, max(tab) - min(tab) + 1)
  s3 <- rep(0, max(tab) - min(tab) + 1)

  # Function to detect unbalanced data
  check_balance <- function(tab) {
    # Check if all elements in `tab` are the same
    length(unique(tab)) == 1
  }


  # Balancing procedure for unbalanced data

  # Calculate items to delete for each group size
  for (i in min(tab):max(tab)) {
    delta <- tab - i
    delta2 <- tab - i
    delta3 <- rep(0,length(delta))
    for (j in 1:length(tab)) {
      if (delta[j] < 0) { # this means we need to delete the whole group
        delta[j] = punish.coeff*tab[j]  # with punishment
        delta2[j] = tab[j] # without punishment
        delta3[j] = 1
      }
    }
    s[i - min(tab) + 1] <- sum(delta) # Total items to delete for this group size with punishment
    s2[i - min(tab) + 1] <- sum(delta2) # without punishment
    s3[i - min(tab) + 1] <- sum(delta3) # total number of group we delete
  }

  # Create a summary data frame for group sizes and items to delete
  S <- data.frame(
    Group_size = c(min(tab):max(tab)),
    Items_to_delete_with_punishment = s,
    Items_to_delete_without_punishment = s2,
    Groups_to_delete= s3
  )

  # View the data frame in R
  # View(S)

  # Display the minimum items to delete, if the data is unbalanced
  # Find the minimum entries to delete with punishment, and for this case show the real number of entries to delete (without punishment)
  if (check_balance(tab)==FALSE){

    min_items_to_delete_without_punishment <- S$Items_to_delete_without_punishment[which.min(S$Items_to_delete_with_punishment)]
    # cat("\nOptimal number of entries to delete for balance:", min_items_to_delete_without_punishment, "\n")

    min_groups_to_delete <- S$Groups_to_delete [which.min(S$Items_to_delete_with_punishment)]
    # cat("\nOptimal number of groups to delete for balance:", min_groups_to_delete, "\n")
  }

  # If data is unbalanced issue a warning with number of entries, groups to delete
  if (check_balance(tab)==FALSE) {
    percentage_items_deleted <- 100*min_items_to_delete_without_punishment / nrow(A)
    
    warning(sprintf(" Your data is unbalanced. Balancing procedure was used and %.1f%% of data was deleted.\n   Deleted entries: %d\n   Deleted groups: %d", percentage_items_deleted, min_items_to_delete_without_punishment, min_groups_to_delete))


    if (balancing.limit != 0.2)
    {
      warning(" You changed the balancing limit to fine-tune the balancing.\n   Increasing the balancing limit might result in the loss of data connections.\n")
    }
  }

  # Check if the ratio of minimum items to delete is within the balancing.limit
  if (check_balance(tab)==FALSE) {
    if (min_items_to_delete_without_punishment / nrow(A) > balancing.limit || min_groups_to_delete / group_num > balancing.limit ) {
      
      balancing.percentage <-balancing.limit*100 # to display balancing limit in percents
      
      stop(sprintf("The share of data that should be deleted is %.1f%% and exceeds the balancing.limit, %s%%.\n\n  Data may not be balanced. If you want to run, adjust the balancing.limit. \n\n", percentage_items_deleted, format(balancing.percentage, trim = TRUE)))

      return(A)  # Optionally return the unbalanced matrix
    }

    # Balance the model matrix A based on optimal size (number of groups and elements in the group)
    target_group_size <- S$Group_size[which.min(S$Items_to_delete_with_punishment)]

    for (group_name in names(group_counts)) {
      group_indices <- which(A$group == group_name)

      if (length(group_indices) > target_group_size) {
        # If the group size is larger than the target, remove excess elements
        to_remove <- length(group_indices) - target_group_size
        group_indices_to_remove <- sample(group_indices, to_remove)
        A <- A[-group_indices_to_remove, ]

      } else if (length(group_indices) < target_group_size) {
        # If the group size is smaller than the target, remove the entire group
        group_indices_to_remove <- group_indices
        A <- A[-group_indices_to_remove, ]
      }
    }
    
      #  delete groups that are unused for the case when groups are given by sequential factor type of the data
      
      if (is.factor(A$group)) {
        A$group <- droplevels(A$group)
      }
    
    
      # Print final counts after balancing
      balanced_group_counts <- table(A$group)
      
      
  #    cat("\nGroups and their sizes after balancing:\n")
  #    print(balanced_group_counts)

    # The balanced model matrix and group sizes for further use
    #return(A)
    group_num <- length(balanced_group_counts)
    group_size <- unique(balanced_group_counts)
    if (length(group_size) > 1){
      warning("Data was not balanced correctly.")
    }
    
  } else {
    group_num <- length(group_counts)
    group_size <- unique(group_counts)
  }
  

  if(FALSE){
   #Use na.omit to remove rows with NA values in relevant columns
    relevant_vars <- c(response_var, predictor_var, control_vars)
    complete_data <- na.omit(data[relevant_vars]) # complete_data <- data[complete.cases(data[relevant_vars]), ]
  
    #Recalculate number of rows
    n_rows <- nrow(complete_data)
  
    #Ensure we have enough data points after removing NAs
    if (n_rows < 1) {
      stop("No complete cases available after removing missing values.")
    }
  
  }

  # Store original data before processing
  original_data <- data
  
  # changing the model matrix into a data frame called data
  data <-as.data.frame(A)

  predictor_var <- names(data)[2]
  
  # check that there are control variables in the data
  if (ncol(data) > 3) {
    control_vars <- names(data)[3:(length(names(data)) - 1)] # do not include last column - it defines groups
  }
  # Extract relevant data columns
  y <- data[[response_var]]
  x <- data[[predictor_var]]
  # if control variables exist
  if (exists("control_vars")) {
    C <- data[control_vars]
    C <- as.matrix(C)
  }
  # Ensure that all variables (response, predictor, and controls) are numeric
  # Check response variable
  if (!is.numeric(y)) {
    stop(paste("The response variable", response_var, "must be numeric."))
  }

  # Check predictor variable
  if (!is.numeric(x)) {
    stop(paste("The predictor variable", predictor_var, "must be numeric."))
  }

  # Check all control variables (if any exist)
  if (exists("control_vars")) {
    nonnum_values <- sapply(C, function(x) !is.numeric(x))

    if (any(nonnum_values)) {
      stop(paste("The following control variables contain non-numeric values:",
                 paste(colnames(C[nonnum_values]), collapse = ", ")))
    }
  }

  if (sum(is.na(y))>0) {
    stop(paste("The response variable", response_var," has null (nan) values."))
  }

  # Check predictor variable
  if (sum(is.na(x))>0) {
    stop(paste("The predictor variable", predictor_var, "has null (nan) values."))
  }
  exists("control_vars")
  # Checking for missing values in control variables
  if (exists("control_vars")) {
    missing_values <- sapply(C, anyNA)

    if (any(missing_values)) {
      stop(paste("The following control variables contain missing values:",
                 paste(colnames(C[missing_values]), collapse = ", ")))
    }
  }

  # If all checks pass, proceed with the function...
  # message("All checks passed. Function is ready to proceed.")


  # Calculate group size (n) and total observations (kn)
  n <- group_size  # Size of each group
  k <- group_num # Number of groups
  kn <- k * n  # Total observations
  
  # Prepare data_CV list
  data_CV <- list(
    y = y,
    x = x,
    k = k,
    n = n,
    kn = kn
  )
  
  # Check if C exists
  if (exists("control_vars")) {
    data_CV$C <- C  # Add C to the list
    data_CV$kc <- length(control_vars)  # Add number of control variables
  }

  
  ML <- estimate_ML_CV(data_CV) # run ML estimator and get a preliminary estimation of b_b for estimate_Bay_CV
  #check if there is any variation between group
  if (ML$tau_x2==0 || ML$tau_yx==0){
    stop(sprintf("The provided data does not include between-group variation for the given model.\n  A two-level model is no longer necessary because the grouping structure does not add information.\n  Consider using a single-level model, treating all observations as independent."))
  }
  #check if there is any variation within group
  if (ML$sigma_x2==0 || ML$sigma_yx==0){
    stop(sprintf("The provided data does not include within-group variation for the given model.\n  All observations within a group are identical (or perfectly correlated).\n  Consider using a simpler model."))
  }
  
  
  data_CV$b_b = ML$beta_b_ML_CV # dummy real value of beta_b
  
  # Call estimate_Bay_CV function with data_CV
  Bay <- estimate_Bay_CV(data_CV)
  
  # recalculate SE of Bayesian estimator with jackknife if a
  if (jackknife == TRUE){
    Bay_jackknife    <- estimate_Bay_CV_SE_jackknife_individual(data_CV)
    Bay$SE_beta_Bay  <- Bay_jackknife$SE_beta_Bay_ML_jackknife_individual
    # Bay$SE_beta_ML   <- Bay_jackknife$SE_beta_ML_jackknife_individual # open in case jackknife for ML needed
    # If there were any gamma-covariates, copy over their SEs
    if (!is.null(Bay_jackknife$SE_gamma_jackknife_individual)) {
      Bay$SE_gamma   <- Bay_jackknife$SE_gamma_jackknife_individual
    }
    
  }

  # Generate the result output

  # Number of control variables (kc)
  if (exists("control_vars")){
    kc <- data_CV$kc
  } else {
    kc<-0
  }
  
  # If control variables are present
  if (kc>0) {
    # Create the list of estimated values dynamically
    Coefficients <- data.frame(
      beta_b = Bay$beta_b_Bay,
      t(sapply(1:kc, function(i) Bay$gamma[i]))   # Dynamic number of gamma columns
    )
    
    colnames(Coefficients) <- c("beta_b", paste0("gamma_", control_vars))  # Adjust gamma column names
  } else {
    Coefficients <- data.frame(
      beta_b = Bay$beta_b_Bay
    )
    colnames(Coefficients) <- c("beta_b")  # Adjust column names
  }
  
  
  if (kc>0) {
    Standard_Error <- data.frame(
      beta_b = Bay$SE_beta_Bay,
      t(sapply(1:kc, function(i) Bay$SE_gamma[i]))
    )
  
    colnames(Standard_Error) <- c("beta_b", paste0("gamma_", control_vars))
  } else {
    Standard_Error <- data.frame(
      beta_b = Bay$SE_beta_Bay
    )
    
    colnames(Standard_Error) <- c("beta_b")
  }
    
  
  if (kc>0) {
    Confidence_Interval <- data.frame(
      Lower = c(Bay$beta_b_Bay - stats::qnorm(1-conf.level/2) * Bay$SE_beta_Bay, Bay$gamma - stats::qnorm(1-conf.level/2) * Bay$SE_gamma),
      Upper = c(Bay$beta_b_Bay + stats::qnorm(1-conf.level/2) * Bay$SE_beta_Bay, Bay$gamma + stats::qnorm(1-conf.level/2) * Bay$SE_gamma)
    )
    
    rownames(Confidence_Interval) <- c("beta_b", paste0("gamma_", control_vars))
  } else {
    Confidence_Interval <- data.frame(
      Lower = c(Bay$beta_b_Bay - stats::qnorm(1-conf.level/2) * Bay$SE_beta_Bay),
      Upper = c(Bay$beta_b_Bay + stats::qnorm(1-conf.level/2) * Bay$SE_beta_Bay)
    )
    
    rownames(Confidence_Interval) <- c("beta_b")
  }

  Confidence_level <- paste0((1 - conf.level) * 100, "%")
  
  if (kc>0) {
    Z_value <- data.frame(
      beta_b = Bay$beta_b_Bay / Bay$SE_beta_Bay,
      t(sapply(1:kc, function(i) Bay$gamma[i] / Bay$SE_gamma[i]))
    )
  
    colnames(Z_value) <- c("beta_b", paste0("gamma_", control_vars))
  } else {
    Z_value <- data.frame(
      beta_b = Bay$beta_b_Bay / Bay$SE_beta_Bay
    )
    
    colnames(Z_value) <- c("beta_b")
  }
  
  if (kc>0) {
    p_value <- data.frame(
      beta_b = 2 * (1 - stats::pnorm(abs(Bay$beta_b_Bay / Bay$SE_beta_Bay))),
      t(sapply(1:kc, function(i) 2 * (1 - stats::pnorm(abs(Bay$gamma[i] / Bay$SE_gamma[i])))))
    )
  
    colnames(p_value) <- c("beta_b", paste0("gamma_", control_vars))
  } else {
    p_value <- data.frame(
      beta_b = 2 * (1 - stats::pnorm(abs(Bay$beta_b_Bay / Bay$SE_beta_Bay)))
    )
    
    colnames(p_value) <- c("beta_b")
  }

  # Create the dynamic call_info string
  call_info <- paste0("mlob(", deparse(formula), ", data = ", data_name, ", group = ", group)

  # Add balancing.limit if it is not default value
  if (!missing(balancing.limit) && balancing.limit != 0.2) {
    call_info <- paste0(call_info, ", balancing.limit = ", balancing.limit)
  }

  # Conditionally add `conf.level` if it's not the default value
  if (!missing(conf.level) && conf.level != 0.05) {
    call_info <- paste0(call_info, ", conf.level = ", 1 - conf.level)
  }

  # Conditionally add `jackknife` if it's not TRUE
  if (!missing(jackknife) && !jackknife) {
    call_info <- paste0(call_info, ", jackknife = ", jackknife)
  }

  call_info <- paste0(call_info, ")")

  # Create call_args list for easier updating
  # Store the original data before it was processed/balanced
  call_args <- list(
    formula = formula,
    data = original_data,  # Store original data, not processed data
    group = group,
    balancing.limit = balancing.limit,
    conf.level = 1 - conf.level,  # Store original confidence level
    jackknife = jackknife,
    punish.coeff = punish.coeff
  )

  result <- list(
    Coefficients = Coefficients,
    Standard_Error = Standard_Error,
    Confidence_Interval = Confidence_Interval,
    Confidence_level = Confidence_level,
    Z_value = Z_value,
    p_value = p_value,
    call_info = call_info,
    call_args = call_args
  )
  
  class(result) <- "mlob_result"  # Assign custom class
  
  # return(result)

  # For unoptimized estimator  estimate_ML
  
  if (kc>0) {
    Coefficients_ML <- data.frame(
      beta_b = Bay$beta_b_ML,  # Using beta_b_ML here
      t(sapply(1:kc, function(i) Bay$gamma[i]))   # Dynamic number of gamma columns
    )
    colnames(Coefficients_ML) <- c("beta_b_ML", paste0("gamma_", control_vars))
  
    Standard_Error_ML <- data.frame(
      beta_b = Bay$SE_beta_ML,
      t(sapply(1:kc, function(i) Bay$SE_gamma[i]))
    )
    colnames(Standard_Error_ML) <- c("beta_b_ML", paste0("gamma_", control_vars))
  
    Confidence_Interval_ML <- data.frame(
      Lower = c(Bay$beta_b_ML - stats::qnorm(1 - conf.level / 2) * Bay$SE_beta_ML, Bay$gamma - stats::qnorm(1 - conf.level / 2) * Bay$SE_gamma),
      Upper = c(Bay$beta_b_ML + stats::qnorm(1 - conf.level / 2) * Bay$SE_beta_ML, Bay$gamma + stats::qnorm(1 - conf.level / 2) * Bay$SE_gamma)
    )
    rownames(Confidence_Interval_ML) <- c("beta_b_ML", paste0("gamma_", control_vars))
  
    Confidence_level_ML <- paste0((1 - conf.level) * 100, "%")
  
    Z_value_ML <- data.frame(
      beta_b = Bay$beta_b_ML / Bay$SE_beta_ML,
      t(sapply(1:kc, function(i) Bay$gamma[i] / Bay$SE_gamma[i]))
    )
    colnames(Z_value_ML) <- c("beta_b_ML", paste0("gamma_", control_vars))
  
    p_value_ML <- data.frame(
      beta_b = 2 * (1 - stats::pnorm(abs(Bay$beta_b_ML / Bay$SE_beta_ML))),
      t(sapply(1:kc, function(i) 2 * (1 - stats::pnorm(abs(Bay$gamma[i] / Bay$SE_gamma[i])))))
    )
    colnames(p_value_ML) <- c("beta_b", paste0("gamma_", control_vars))
    
  } else {
    
    Coefficients_ML <- data.frame(
      beta_b = Bay$beta_b_ML  # Using beta_b_ML here
    )
    
    colnames(Coefficients_ML) <- c("beta_b_ML")
    
    Standard_Error_ML <- data.frame(
      beta_b = Bay$SE_beta_ML
    )
    
    colnames(Standard_Error_ML) <- c("beta_b_ML")
    
    Confidence_Interval_ML <- data.frame(
      Lower = c(Bay$beta_b_ML - stats::qnorm(1 - conf.level / 2) * Bay$SE_beta_ML),
      Upper = c(Bay$beta_b_ML + stats::qnorm(1 - conf.level / 2) * Bay$SE_beta_ML)
    )
    
    rownames(Confidence_Interval_ML) <- c("beta_b_ML")
    
    Confidence_level_ML <- paste0((1 - conf.level) * 100, "%")
    
    Z_value_ML <- data.frame(
      beta_b = Bay$beta_b_ML / Bay$SE_beta_ML
    )
    
    colnames(Z_value_ML) <- c("beta_b_ML")
    
    p_value_ML <- data.frame(
      beta_b = 2 * (1 - stats::pnorm(abs(Bay$beta_b_ML / Bay$SE_beta_ML)))
    )
    
    colnames(p_value_ML) <- c("beta_b")
  }

  result$Coefficients_ML = Coefficients_ML
  result$Standard_Error_ML = Standard_Error_ML
  result$Confidence_Interval_ML = Confidence_Interval_ML
  result$Confidence_level_ML = Confidence_level_ML
  result$Z_value_ML = Z_value_ML
  result$p_value_ML = p_value_ML

  return(result)

}