#' Print Method for Design Output from ForLion Algorithms
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

  ####identify the variable names
  num_vars <- ncol(x_factor)
  if (!is.null(x$var.names) && length(x$var.names) == num_vars) {
    var_names <- x$var.names
  } else {
    var_names <- paste0("X", 1:num_vars) # default the x1,x2,x3...
  }

  # Determine column widths dynamically
  # max_widths <- apply(x_factor, 2, function(col) max(nchar(sprintf("%.4f", col)), na.rm = TRUE))
  # max_width_total <- sum(max_widths) + length(max_widths) - 1  # Account for commas
  # max_alloc_width <- max(nchar(sprintf("%.4f", p)))
  data_widths <- apply(x_factor, 2, function(col) max(nchar(sprintf("%.4f", col)), na.rm = TRUE))
  header_widths <- nchar(var_names)
  final_var_widths <- pmax(data_widths, header_widths)

  count_width <- max(nchar("Count"), nchar(nrow(x_factor)))
  alloc_width <- max(nchar("Allocation"), nchar(sprintf("%.4f", p)))


  cat("Design Output\n")
  #cat(paste0(rep("=", 6 + 2 + max_width_total + 2 + max_alloc_width+8), collapse = ""), "\n")
  total_width <- count_width + sum(final_var_widths) + alloc_width + (num_vars + 1) * 2
  cat(paste0(rep("=", total_width), collapse = ""), "\n")

  # Print table header
  cat(sprintf("%-*s", count_width, "Count")) # print count
  cat("  ")

  for (i in seq_along(var_names)) {
    cat(sprintf("%-*s", final_var_widths[i], var_names[i]))
    cat("  ") #
  }

  cat(sprintf("%s\n", "Allocation"))


  # cat(sprintf("%-6s  %-*s  %-*s\n",
  #             "Count", max_width_total, "Design Point(s)", max_alloc_width, "Allocation"))
  cat(paste0(rep("-", total_width), collapse = ""), "\n")


  # Print the table rows with aligned columns
  for (i in seq_len(nrow(x_factor))) {
    # design_point <- paste(sprintf(paste0("%", max_widths, ".4f"), x_factor[i, ]), collapse = ", ")  # Align design points
    # allocation <- sprintf("%.*f", 4, p[i])                                     # Round allocation
    # cat(sprintf("%-6d  %-*s  %-*s\n", i, max_width_total, design_point, max_alloc_width, allocation))

    cat(sprintf("%-*d", count_width, i))
    cat("  ")

    for (j in seq_along(final_var_widths)) {
      number_as_string <- sprintf("%.4f", x_factor[i, j])
      cat(sprintf("%*s", final_var_widths[j], number_as_string))
      cat("  ")
    }

    cat(sprintf("%.4f\n", p[i]))
  }

  cat(paste0(rep("=", total_width), collapse = ""), "\n")

  # Print other elements in the list
  other_elements <- setdiff(names(x), c("x.factor", "p", "var.names"))
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
