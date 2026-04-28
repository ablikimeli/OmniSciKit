#' Convert Wide to Long Format
#'
#' Convert data from wide to long format
#' @param data Data frame in wide format
#' @param id_cols Columns to keep as identifiers
#' @param names_to Name of new column for variable names
#' @param values_to Name of new column for values
#' @return Data frame in long format
#' @examples
#' data(omni_depression)
#' long_data <- convert_wide_to_long(
#'   omni_depression,
#'   id_cols = c("patient_id", "group"),
#'   names_to = "timepoint",
#'   values_to = "score"
#' )
#' @export
convert_wide_to_long <- function(data, id_cols, names_to = "variable", values_to = "value") {
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    # Fallback using base R
    id_data <- data[, id_cols, drop = FALSE]
    measure_cols <- setdiff(names(data), id_cols)

    result <- data.frame()
    for (col in measure_cols) {
      temp <- id_data
      temp[[names_to]] <- col
      temp[[values_to]] <- data[[col]]
      result <- rbind(result, temp)
    }

    return(result)
  }

  tidyr::pivot_longer(
    data,
    cols = -dplyr::all_of(id_cols),
    names_to = names_to,
    values_to = values_to
  )
}

#' Convert Long to Wide Format
#'
#' Convert data from long to wide format
#' @param data Data frame in long format
#' @param id_cols Columns to keep as identifiers
#' @param names_from Column containing variable names
#' @param values_from Column containing values
#' @return Data frame in wide format
#' @examples
#' # First create long data, then convert back
#' data(omni_depression)
#' @export
convert_long_to_wide <- function(data, id_cols, names_from = "variable", values_from = "value") {
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    # Fallback using base R reshape
    return(stats::reshape(
      data,
      idvar = id_cols,
      timevar = names_from,
      direction = "wide",
      v.names = values_from
    ))
  }

  tidyr::pivot_wider(
    data,
    id_cols = dplyr::all_of(id_cols),
    names_from = dplyr::all_of(names_from),
    values_from = dplyr::all_of(values_from)
  )
}

#' Standardize Variables
#'
#' Standardize numeric variables to mean=0, sd=1
#' @param data Data frame
#' @param vars Variables to standardize (NULL = all numeric)
#' @return Data frame with standardized variables
#' @examples
#' data(omni_glucose)
#' std_data <- standardize_vars(omni_glucose, vars = c("age", "bmi"))
#' @export
standardize_vars <- function(data, vars = NULL) {
  if (is.null(vars)) {
    vars <- names(data)[sapply(data, is.numeric)]
  }

  result <- data
  for (var in vars) {
    if (is.numeric(data[[var]])) {
      result[[var]] <- scale(data[[var]])[, 1]
    }
  }

  return(result)
}

#' Merge Datasets
#'
#' Merge two data frames by common columns
#' @param x First data frame
#' @param y Second data frame
#' @param by Columns to merge by
#' @param all Logical for full outer join
#' @return Merged data frame
#' @examples
#' data(omni_glucose)
#' # Merge example (same data for demonstration)
#' merged <- merge_datasets(omni_glucose[1:10, ], omni_glucose[11:20, ], by = "patient_id")
#' @export
merge_datasets <- function(x, y, by = NULL, all = FALSE) {
  if (is.null(by)) {
    by <- intersect(names(x), names(y))
    if (length(by) == 0) {
      stop("No common columns found for merging")
    }
  }

  merge(x, y, by = by, all = all)
}

#' Import CSV File
#'
#' Import data from CSV file with automatic encoding detection
#' @param file Path to CSV file
#' @param header Logical for header row
#' @param sep Separator character
#' @param encoding File encoding
#' @return Data frame
#' @examples
#' # import_csv("data.csv")
#' @export
import_csv <- function(file, header = TRUE, sep = ",", encoding = "UTF-8") {
  utils::read.csv(file, header = header, sep = sep, fileEncoding = encoding,
                   stringsAsFactors = FALSE)
}

#' Export CSV File
#'
#' Export data frame to CSV file
#' @param data Data frame
#' @param file Output file path
#' @param row.names Include row names
#' @param encoding File encoding
#' @return NULL
#' @examples
#' # export_csv(omni_glucose, "output.csv")
#' @export
export_csv <- function(data, file, row.names = FALSE, encoding = "UTF-8") {
  utils::write.csv(data, file, row.names = row.names, fileEncoding = encoding)
}
