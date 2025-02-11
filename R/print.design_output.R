#' Print Method for Design Output from ForLion Algorithm
#'
#' @description Custom print method for a list containing design information.
#' @param x An object of class `design_output`.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns `x`.
#' @export
#' @method print design_output
#'
#'

print.design_output <- function(x, ...) {
  if (!inherits(x, "design_output")) stop("Object is not of class 'design_output'.")

  # Validate and process x.factor and p
  x_factor <- x$x.factor
  p <- x$p

  # Handle vector and matrix cases for x.factor
  if (is.vector(x_factor)) {
    x_factor <- matrix(x_factor, ncol = 1)  # Convert vector to single-column matrix
  }

  if (!is.matrix(x_factor) || length(p) != nrow(x_factor)) {
    stop("Invalid input: $x.factor must be a matrix, and $p must have the same length as rows in $x.factor.")
  }

  # Determine column widths dynamically
  max_widths <- apply(x_factor, 2, function(col) max(nchar(sprintf("%.4f", col)), na.rm = TRUE))
  max_width_total <- sum(max_widths) + length(max_widths) - 1  # Account for commas
  max_alloc_width <- max(nchar(sprintf("%.4f", p)))

  cat("Design Output\n")
  cat(paste0(rep("=", 6 + 2 + max_width_total + 2 + max_alloc_width+8), collapse = ""), "\n")

  # Print table header
  cat(sprintf("%-6s  %-*s  %-*s\n",
              "Count", max_width_total, "Design Point(s)", max_alloc_width, "Allocation"))
  cat(paste0(rep("-", 6 + 2 + max_width_total + 2 + max_alloc_width+8), collapse = ""), "\n")


  # Print the table rows with aligned columns
  for (i in seq_len(nrow(x_factor))) {
    design_point <- paste(sprintf(paste0("%", max_widths, ".4f"), x_factor[i, ]), collapse = ", ")  # Align design points
    allocation <- sprintf("%.*f", 4, p[i])                                     # Round allocation
    cat(sprintf("%-6d  %-*s  %-*s\n", i, max_width_total, design_point, max_alloc_width, allocation))
  }

  cat(paste0(rep("=", 6 + 2 + max_width_total + 2 + max_alloc_width+8), collapse = ""), "\n")

  # Print other elements in the list
  other_elements <- setdiff(names(x), c("x.factor", "p"))
  for (name in other_elements) {
    cat("\n", name, ":\n", sep = "")
    element <- x[[name]]
    det_name<-"det"
    rel.effi_name<-"rel.efficiency"
    if(name==det_name || name==rel.effi_name){
      element=element
    }else if (is.numeric(element)) {
      element <- round(element, 4)  # Round numeric values to 4 decimal places
    }
    print(element)
  }

  invisible(x)
}







