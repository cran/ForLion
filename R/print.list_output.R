#' Print Method for list_output Objects
#'
#' @description Custom print method for objects of class `list_output`.
#' @param x An object of class `list_output`.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns `x` (the input object).
#' @export
#' @method print list_output
#'
#'

print.list_output <- function(x, ...) {
  if (!is.list(x)) {
    stop("The object is not a list.")
  }

  cat("D-Optimal Design Results:\n")
  cat(rep("=", 80), "\n", sep = "")

  for (i in seq_along(x)) {
    element <- x[[i]]
    element_name <- names(x)[i]

    cat(element_name, ":\n")

    if (is.matrix(element)) {
      # Pretty print for matrices
      print(format(element, justify = "right"), quote = FALSE)
    } else if (is.data.frame(element)) {
      # Pretty print for data frames
      print(element, row.names = FALSE)
    } else if (is.numeric(element)) {
      # Format numeric values dynamically
      formatted <- sapply(element, function(value) {
        if (value == round(value, digits = 0) && grepl("\\.0$", format(value, nsmall = 1))) {
          # If the value is an integer with a ".0", keep it as is
          format(value, nsmall = 1)
        } else if (value == round(value)) {
          # If the value is a true integer, print it as an integer
          as.character(value)
        } else if (abs(value) < 1e-3 || abs(value) >= 1e+3) {
          # For very small or large values, use scientific notation with 4 decimal places
          formatC(value, format = "e", digits = 4)
        } else {
          # For other values, show up to 4 decimal places, removing unnecessary trailing zeros
          sub("\\.?0+$", "", formatC(value, format = "f", digits = 4))
        }
      })
      cat(paste(formatted, collapse = ", "), "\n")
    } else if (is.logical(element)) {
      # Print logical vectors
      cat(paste(element, collapse = ", "), "\n")
    } else if (is.character(element)) {
      # Print character vectors
      cat(paste0('"', element, '"', collapse = ", "), "\n")
    } else if (is.list(element)) {
      # Indicate nested lists
      cat("Nested list with", length(element), "elements\n")
    } else {
      # Fallback for unknown types
      cat("Unsupported element type:", class(element), "\n")
    }

    cat(rep("-", 80), "\n", sep = "")
  }

  invisible(x)
}
