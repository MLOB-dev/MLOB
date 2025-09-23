#' S3 Methods for mlob_result Class
#'
#' This file contains all S3 methods for the mlob_result class
#' from the MultiLevelOptimalBayes package.
#'
#' @keywords internal

#' Print method for mlob_result objects
#'
#' @param x An object of class mlob_result
#' @param ... Additional arguments passed to print
#' @export
print.mlob_result <- function(x, ...) {
  # Custom print method for mlob_result
  # Extract call information
  cat("Call:\n", x$call_info, "\n\n")

  # Print each section without row names
  cat("Coefficients\n")
  print(x$Coefficients, row.names = FALSE)

  cat("\nStandard_Error\n")
  print(x$Standard_Error, row.names = FALSE)

  cat("\nConfidence_Interval (", x$Confidence_level, ")\n", sep = "")
  print(x$Confidence_Interval)

  cat("\nZ value\n")
  print(x$Z_value, row.names = FALSE)

  cat("\np value\n")
  print(x$p_value, row.names = FALSE)
  
  invisible(x)
}

#' Summary method for mlob_result objects
#'
#' @param object An object of class mlob_result
#' @param ... Additional arguments passed to summary
#' @export
summary.mlob_result <- function(object, ...) {
  # Custom print method for summary_mlob_result
  # Function to assign stars based on p-values
  signif_stars <- function(pval) {
    if (pval < 0.001) {
      return("***")
    } else if (pval < 0.01) {
      return("**")
    } else if (pval < 0.05) {
      return("*")
    } else if (pval < 0.1) {
      return(".")
    } else {
      return(" ")
    }
  }


  # Prepare the summary table with significance stars
  summary_table <- data.frame(
    Estimate = as.numeric(object$Coefficients),
    `Std. Error` = as.numeric(object$Standard_Error),
    `Lower CI` = as.numeric(object$Confidence_Interval$Lower),
    `Upper CI` = as.numeric(object$Confidence_Interval$Upper),
    `Z value` = as.numeric(object$Z_value),
    `Pr(>|z|)` = as.numeric(object$p_value),
    Significance = sapply(object$p_value, signif_stars)
  )

  # If control variables are present in the model
  if (ncol(object$Coefficients)>1){
  # Extract control variable names from the Coefficients (removing the 'beta_b' column)
  control_vars <- colnames(object$Coefficients)[-1]  # Get control variable names (gamma_Age, gamma_Education, etc.)

  # Update the row names for the summary table using the actual control variable names
  rownames(summary_table) <- c("beta_b", control_vars)
  } else {
    # Update the row names for the summary table using the actual control variable names
    rownames(summary_table) <- c("beta_b")
  }
  
  # Redact column names for consistency
  colnames(summary_table) <- c("Estimate", "Std. Error", paste0("Lower CI (", object$Confidence_level, ")"),
                               paste0("Upper CI (", object$Confidence_level, ")"), "Z value", "Pr(>|z|)", "Significance")


  # Print the summary header
  cat("Call:\n", object$call_info, "\n\n")
  cat("Summary of Coefficients:\n")

  print(summary_table, row.names = TRUE)


  invisible(summary_table)  # Return the table invisibly so it doesn't clutter the output


  # Prepare the summary table For unoptimized estimator ML with significance stars
  summary_table_ML<- data.frame(
    Estimate = as.numeric(object$Coefficients_ML),
    `Std. Error` = as.numeric(object$Standard_Error_ML),
    `Lower CI` = as.numeric(object$Confidence_Interval_ML$Lower),
    `Upper CI` = as.numeric(object$Confidence_Interval_ML$Upper),
    `Z value` = as.numeric(object$Z_value_ML),
    `Pr(>|z|)` = as.numeric(object$p_value_ML),
    Significance = sapply(object$p_value_ML, signif_stars)
  )

  ## Extract control variable names from the Coefficients (removing the 'beta_b' column)
  # control_vars <- colnames(object$Coefficients)[-1]  # Get control variable names (gamma_Age, gamma_Education, etc.)
  
  # If control variables are present in the model
  if (ncol(object$Coefficients)>1){
    # Define the number of rows based on the number of coefficients (beta_b + gamma terms)
    rownames(summary_table_ML) <- c("beta_b", control_vars)
  } else {
    rownames(summary_table_ML) <- c("beta_b")
  }
  
  # Redact column names for consistency
  colnames(summary_table_ML) <- c("Estimate", "Std. Error", paste0("Lower CI (", object$Confidence_level, ")"),
                               paste0("Upper CI (", object$Confidence_level, ")"), "Z value", "Pr(>|z|)", "Significance")


  # Print the ML part of summary

  cat("\n\nFor comparison, summary of coefficients from unoptimized analysis (ML):\n")

  print(summary_table_ML, row.names = TRUE)

  # Add significance codes
  cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  
  
  # Extract first SE of ML and of regularized Bayes (or NA if missing)
  se_ml  <- if (!is.null(object$Standard_Error_ML) &&
                length(object$Standard_Error_ML) >= 1) object$Standard_Error_ML[[1]] else NA_real_
  se_reg <- if (!is.null(object$Standard_Error) &&
                length(object$Standard_Error) >= 1) object$Standard_Error[[1]] else NA_real_
  
  # Compute ratio (will be NaN/Inf if se_reg is zero or NA)
  se_ratio <- se_ml / se_reg
  
  # Display the note header
  cat("\nNote:\n")
  
  # Three-way branching without any stop()
  if (!is.finite(se_ratio)) {
    
    # Case 1: missing, infinite, zero, etc.
    cat(
      "  The standard errors cannot be compared because one or both are missing, zero, NA, NaN, or infinite.\n",
      sep = ""
    )
    
  } else if (se_ratio >= 0.9 & se_ratio <= 1.1) {
    
    # Case 2: roughly equal
    cat(
      "  The standard error from unoptimized ML estimation is approximately the same as the standard error obtained through our optimization procedure,\n",
      "  meaning that both approaches yield similarly accurate estimates.\n",
      sep = ""
    )
    
  } else if (se_ratio > 1.1) {
    
    pct_diff <- (se_ratio - 1) * 100
    
    # pick digits: none if large, 1 if moderate, 2 if small
    digits <- if (pct_diff >= 200) {
      0
    } else if (pct_diff >= 100) {
      1
    } else if (pct_diff >= 50) {
      1
    } else {
      2
    }
    
    # format with exactly that many decimals
    fmt_pct <- format(round(pct_diff, digits), nsmall = digits)
    
    cat(
      sprintf(
        "  The standard error from unoptimized ML estimation is about %s%% larger than the standard error obtained through our optimization procedure,\n",
        fmt_pct
      ),
      "  meaning that the optimized estimates are more accurate.\n",
      "  Concerning the estimates themselves, the unoptimized ML estimates may\n",
      "  differ greatly from the optimized estimates and should not be reported.\n",
      "  As the optimized estimates are always at least as accurate as the\n",
      "  unoptimized ML estimates,\n",
      "  please use them and their corresponding standard errors (first table of\n",
      "  output) for interpretation and reporting.\n",
      "  For more information, see Dashuk et al. (2025a).\n",
      sep = ""
    )
  }
}

#' Extract coefficients from mlob_result objects
#'
#' @param object An object of class mlob_result
#' @param ... Additional arguments passed to coef
#' @return A data frame with coefficients
#' @export
coef.mlob_result <- function(object, ...) {
  return(object$Coefficients)
}

#' Generic function to extract standard errors
#'
#' @param object An object
#' @param ... Additional arguments
#' @export
se <- function(object, ...) {
  UseMethod("se")
}

#' Extract standard errors from mlob_result objects
#'
#' @param object An object of class mlob_result
#' @param ... Additional arguments passed to se
#' @return A named vector of standard errors
#' @export
se.mlob_result <- function(object, ...) {
  # Extract standard errors from the Standard_Error component
  if (is.null(object$Standard_Error)) {
    stop("Standard errors not available")
  }
  
  # Extract standard errors
  se_values <- as.numeric(object$Standard_Error[1, ])
  
  # The column names from Standard_Error are already the parameter names
  names(se_values) <- colnames(object$Standard_Error)
  
  return(se_values)
}

#' Extract variances from mlob_result objects
#'
#' @param object An object of class mlob_result
#' @param ... Additional arguments passed to vcov
#' @return A named vector of variances (squared standard errors)
#' @export
vcov.mlob_result <- function(object, ...) {
  # Extract standard errors and convert to variances
  se_values <- se.mlob_result(object)
  
  # Return variances (squared standard errors)
  return(se_values^2)
}

#' Extract confidence intervals from mlob_result objects
#'
#' @param object An object of class mlob_result
#' @param parm Parameters for which to extract confidence intervals
#' @param level Confidence level (can be different from the model's default)
#' @param ... Additional arguments passed to confint
#' @return A matrix with confidence intervals
#' @export
confint.mlob_result <- function(object, parm, level = 0.95, ...) {
  if (missing(parm)) {
    parm <- colnames(object$Coefficients)
  }
  
  # Find the indices of the requested parameters
  param_indices <- match(parm, colnames(object$Coefficients))
  
  if (any(is.na(param_indices))) {
    stop("Some parameters not found in the model")
  }
  
  # Set default level if not provided
  if (is.null(level)) {
    # Convert the model's confidence level from percentage string to numeric
    level <- if (is.character(object$Confidence_level)) {
      as.numeric(gsub("%", "", object$Confidence_level)) / 100
    } else {
      object$Confidence_level
    }
  }
  
  # Validate confidence level
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
    stop("The 'level' argument must be a numeric scalar with value between 0 and 1 (exclusive).")
  }
  
  # Convert the model's confidence level from percentage string to numeric
  model_level <- if (is.character(object$Confidence_level)) {
    as.numeric(gsub("%", "", object$Confidence_level)) / 100
  } else {
    object$Confidence_level
  }
  
  # Check if user requested the default confidence level
  if (abs(level - model_level) < 1e-10) {
    # Return pre-calculated confidence intervals for the requested parameters
    if (is.null(object$Confidence_Interval)) {
      stop("Confidence intervals not available")
    }
    
      # Extract the requested parameters from the pre-calculated intervals
      ci_result <- object$Confidence_Interval[param_indices, , drop = FALSE]
      
      # Set row names to parameter names
      rownames(ci_result) <- parm
      
      # Calculate the actual confidence level percentages for column names
      model_level <- if (is.character(object$Confidence_level)) {
        as.numeric(gsub("%", "", object$Confidence_level)) / 100
      } else {
        object$Confidence_level
      }
      alpha <- 1 - model_level
      lower_percent <- (alpha/2) * 100
      upper_percent <- (1 - alpha/2) * 100
      
      # Update column names to show actual percentages
      colnames(ci_result) <- c(paste0(lower_percent, "%"), paste0(upper_percent, "%"))
      
      return(ci_result)
  }
  
  # If different confidence level requested, recalculate
  # Extract coefficients and standard errors
  if (is.null(object$Coefficients) || is.null(object$Standard_Error)) {
    stop("Coefficients or standard errors not available")
  }
  
  # Get the requested parameters
  coef_values <- as.numeric(object$Coefficients[1, param_indices])
  se_values <- as.numeric(object$Standard_Error[1, param_indices])
  
  # Calculate critical value for the specified confidence level
  # Following the same logic as in mlob_working.R (lines 507-521)
  critical_value <- stats::qnorm((1 - level)/2)
  
  # Calculate confidence intervals: estimate ± critical_value * standard_error
  lower_ci <- coef_values + critical_value * se_values
  upper_ci <- coef_values - critical_value * se_values
  
  # Create the confidence interval matrix
  ci_matrix <- cbind(Lower = lower_ci, Upper = upper_ci)
  rownames(ci_matrix) <- parm
  
  # Calculate the actual confidence level percentages for column names
  # level is the confidence level (e.g., 0.95), so alpha = 1 - level
  alpha <- 1 - level
  lower_percent <- (alpha/2) * 100
  upper_percent <- (1 - alpha/2) * 100
  
  colnames(ci_matrix) <- c(paste0(lower_percent, "%"), paste0(upper_percent, "%"))
  
  return(ci_matrix)
}

#' Convert mlob_result to data frame
#'
#' @param x An object of class mlob_result
#' @param ... Additional arguments passed to as.data.frame
#' @return A data frame representation of the results
#' @export
as.data.frame.mlob_result <- function(x, ...) {
  # Check if required components are available
  if (is.null(x$Coefficients) || is.null(x$Standard_Error)) {
    stop("Coefficients or standard errors not available")
  }
  
  # Create a comprehensive data frame with all results
  result_df <- data.frame(
    Estimate = as.numeric(x$Coefficients[1, ]),
    `Std. Error` = as.numeric(x$Standard_Error[1, ]),
    `Lower CI` = as.numeric(x$Confidence_Interval$Lower),
    `Upper CI` = as.numeric(x$Confidence_Interval$Upper),
    `Z value` = as.numeric(x$Z_value[1, ]),
    `Pr(>|z|)` = as.numeric(x$p_value[1, ]),
    stringsAsFactors = FALSE
  )
  
  # Set row names to parameter names
  rownames(result_df) <- colnames(x$Coefficients)
  
  # Update column names to include confidence level
  colnames(result_df) <- c("Estimate", "Std. Error",
                           paste0("Lower CI (", x$Confidence_level, ")"),
                           paste0("Upper CI (", x$Confidence_level, ")"),
                           "Z value", "Pr(>|z|)")
  
  return(result_df)
}

#' Get dimensions of mlob_result objects
#'
#' @param x An object of class mlob_result
#' @return A vector with dimensions
#' @export
dim.mlob_result <- function(x) {
  return(c(nrow(x$Coefficients), ncol(x$Coefficients)))
}

#' Get length of mlob_result objects
#'
#' @param x An object of class mlob_result
#' @return The number of parameters/coefficients
#' @export
length.mlob_result <- function(x) {
  return(ncol(x$Coefficients))
}

#' Get names of mlob_result objects
#'
#' @param x An object of class mlob_result
#' @return The names of the parameters/coefficients
#' @export
names.mlob_result <- function(x) {
  return(colnames(x$Coefficients))
}

#' Update mlob_result objects
#'
#' @param object An object of class mlob_result
#' @param formula Updated formula
#' @param data Updated data
#' @param group Updated group variable
#' @param balancing.limit Updated balancing limit (default: 0.2)
#' @param conf.level Updated confidence level (default: 0.95)
#' @param jackknife Updated jackknife setting (default: FALSE)
#' @param punish.coeff Updated punishment coefficient (default: 2)
#' @param ... Additional arguments passed to mlob
#' @return Updated mlob_result object
#' @export
update.mlob_result <- function(object, formula = NULL, data = NULL, 
                              group = NULL, balancing.limit = NULL, 
                              conf.level = NULL, jackknife = NULL, 
                              punish.coeff = NULL, ...) {
  # Check if call_args is available
  if (is.null(object$call_args)) {
    stop("Cannot update: call_args not available in the object")
  }
  
  # Get the additional arguments
  dots <- list(...)
  
  # Check if only conf.level is being updated (optimization case)
  # This is the only parameter that can be optimized without full recalculation
  only_conf_level_update <- !is.null(conf.level) && 
                           is.null(formula) && 
                           is.null(data) && 
                           is.null(group) &&
                           is.null(balancing.limit) &&
                           is.null(jackknife) &&
                           is.null(punish.coeff) &&
                           length(dots) == 0
  
  if (only_conf_level_update) {
    # Optimized case: only confidence level changed
    new_conf_level <- conf.level
    
    # Create a copy of the original object
    result <- object
    
    # Update the call_args to reflect the new confidence level
    result$call_args$conf.level <- new_conf_level
    
    # Update the call_info string using the actual formula and parameters
    call_info <- paste0("mlob(", deparse(result$call_args$formula), ", data = data, group = ", result$call_args$group, ", conf.level = ", new_conf_level)
    
    # Add other non-default arguments if they exist
    if (!is.null(result$call_args$balancing.limit) && result$call_args$balancing.limit != 0.2) {
      call_info <- paste0(call_info, ", balancing.limit = ", result$call_args$balancing.limit)
    }
    if (!is.null(result$call_args$jackknife) && !result$call_args$jackknife) {
      call_info <- paste0(call_info, ", jackknife = ", result$call_args$jackknife)
    }
    if (!is.null(result$call_args$punish.coeff) && result$call_args$punish.coeff != 2) {
      call_info <- paste0(call_info, ", punish.coeff = ", result$call_args$punish.coeff)
    }
    
    call_info <- paste0(call_info, ")")
    result$call_info <- call_info
    
    # Update confidence level display
    result$Confidence_level <- paste0((new_conf_level) * 100, "%")
    
    # Recalculate confidence intervals using confint() method
    new_conf_intervals <- stats::confint(object, level = new_conf_level)
    # Convert matrix back to data frame with Lower/Upper columns to maintain compatibility
    result$Confidence_Interval <- data.frame(
      Lower = as.numeric(new_conf_intervals[, 1]),
      Upper = as.numeric(new_conf_intervals[, 2])
    )
    
    # Update ML confidence intervals as well
    if (!is.null(object$Confidence_Interval_ML)) {
      # For ML, we need to recalculate based on the new confidence level
      # Convert confidence level to significance level for internal calculations
      alpha <- 1 - new_conf_level
      
      # Recalculate ML confidence intervals
      if (ncol(object$Coefficients_ML) > 1) {
        # Multiple parameters case
        result$Confidence_Interval_ML <- data.frame(
          Lower = c(object$Coefficients_ML$beta_b_ML - stats::qnorm(1 - alpha/2) * object$Standard_Error_ML$beta_b_ML,
                   sapply(2:ncol(object$Coefficients_ML), function(i) {
                     object$Coefficients_ML[1, i] - stats::qnorm(1 - alpha/2) * object$Standard_Error_ML[1, i]
                   })),
          Upper = c(object$Coefficients_ML$beta_b_ML + stats::qnorm(1 - alpha/2) * object$Standard_Error_ML$beta_b_ML,
                   sapply(2:ncol(object$Coefficients_ML), function(i) {
                     object$Coefficients_ML[1, i] + stats::qnorm(1 - alpha/2) * object$Standard_Error_ML[1, i]
                   }))
        )
        rownames(result$Confidence_Interval_ML) <- c("beta_b_ML", 
                                                    colnames(object$Coefficients_ML)[-1])
      } else {
        # Single parameter case
        result$Confidence_Interval_ML <- data.frame(
          Lower = object$Coefficients_ML$beta_b_ML - stats::qnorm(1 - alpha/2) * object$Standard_Error_ML$beta_b_ML,
          Upper = object$Coefficients_ML$beta_b_ML + stats::qnorm(1 - alpha/2) * object$Standard_Error_ML$beta_b_ML
        )
        rownames(result$Confidence_Interval_ML) <- c("beta_b_ML")
      }
      
      result$Confidence_level_ML <- paste0((1 - alpha) * 100, "%")
    }
    
    return(result)
  }
  
  # Full recalculation case: any parameter other than just conf.level is changed
  # Start with the original call arguments
  call_args <- object$call_args
  
  # Update parameters if provided
  if (!is.null(formula)) {
    call_args$formula <- formula
  }
  if (!is.null(data)) {
    call_args$data <- data
  }
  if (!is.null(group)) {
    call_args$group <- group
  }
  if (!is.null(balancing.limit)) {
    call_args$balancing.limit <- balancing.limit
  }
  if (!is.null(conf.level)) {
    call_args$conf.level <- conf.level
  }
  if (!is.null(jackknife)) {
    call_args$jackknife <- jackknife
  }
  if (!is.null(punish.coeff)) {
    call_args$punish.coeff <- punish.coeff
  }
  
  # Update other arguments from ...
  for (name in names(dots)) {
    call_args[[name]] <- dots[[name]]
  }
  
  # Execute the updated call - ensure data is properly handled
  # Create a proper call structure that mlob() can handle
  # We need to create a call that looks like mlob(formula, data = data_name, group = group, ...)
  
  # Create a temporary environment with the data available
  temp_env <- new.env()
  temp_env$data <- call_args$data
  
  # Build the call as a string to avoid deparsing issues
  call_str <- "mlob(formula = call_args$formula, data = data, group = call_args$group"
  
  # Add optional arguments if they're not default values
  if (!is.null(call_args$balancing.limit) && call_args$balancing.limit != 0.2) {
    call_str <- paste0(call_str, ", balancing.limit = ", call_args$balancing.limit)
  }
  if (!is.null(call_args$conf.level) && call_args$conf.level != 0.95) {
    call_str <- paste0(call_str, ", conf.level = ", call_args$conf.level)
  }
  if (!is.null(call_args$jackknife) && call_args$jackknife != FALSE) {
    call_str <- paste0(call_str, ", jackknife = ", call_args$jackknife)
  }
  if (!is.null(call_args$punish.coeff) && call_args$punish.coeff != 2) {
    call_str <- paste0(call_str, ", punish.coeff = ", call_args$punish.coeff)
  }
  
  call_str <- paste0(call_str, ")")
  
  # Execute the call in the temporary environment
  result <- eval(parse(text = call_str), envir = temp_env)
  
  return(result)
}
