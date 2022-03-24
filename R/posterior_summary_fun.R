# Getting summary statistics
#' @title Calculate summary statistics from Samples
#' @description Takes samples and calculates summary statistics (mean, standard
#' deviation, and quantiles (if desired)).
#' @param .data \code{data.frame} with samples to be summarised.
#' @param probs Percentiles to provide quantiles at. Set to NULL to skip
#' computing quantiles.
#' @importFrom dplyr %>%
#' @export
posterior_summary_fun <- function(.data, probs = c(0.025, 0.5, 0.975)) {
  
  #global bindings for data.table non-standard evaluation
  . <- value <- NULL
  
  probs <- sort(probs)

  # ensure numeric columns are after categorical
  .data <- .data %>%
    dplyr::relocate(
      .data$population | dplyr::contains("samp_"),
      .after = dplyr::everything()
    )

  # pull locations of columns to "group by"
  id_cols <- seq_along(names(.data)[!grepl("samp", names(.data))])

  # use data.table as this can be quite slow for larger countries
  if (!inherits(.data, "data.table")) .data <- data.table::setDT(.data)

  # pivot to calculate row means and sds of samples for each stratification
  .data_long <- data.table::melt(.data,
    id.vars = id_cols,
    measure.vars = c(paste0("samp_", 1:100))
  )
  .data <- .data_long[,
    "."
    (mean = mean(value, na.rm = TRUE),
      sd = stats::sd(value, na.rm = TRUE)),
    keyby = c(names(.data)[id_cols])
  ] # group by all categories]

  # calculate median and CI
  if (!is.null(probs)) {
    quantiles <- .data_long[, {
        quantiles <- stats::quantile(value,
          probs,
          na.rm = TRUE,
          names = FALSE
        )
        ":="
        list(
          lower  = quantiles[1],
          median = quantiles[2],
          upper  = quantiles[3]
        )},
      keyby = c(names(.data)[id_cols]),
    ]

    return(as.data.frame(merge(.data, quantiles)))
  } else {
    return(.data)
  }
}
