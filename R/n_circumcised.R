#' @title Calculate number of people circumcised
#' @description Calculate number of people circumcised (as well as unmet need).
#' @param results Results with samples for number of circumcisions performed
#' in each region.
#' @importFrom dplyr %>%
#' @export 
n_circumcised <- function(results) {
    # Getting number of circumcised men
    tmp <- split(results, results$type)

    # get circumcised population by type
    tmp <- lapply(tmp, function(x) {
        x %>%
            dplyr::mutate(
                dplyr::across(dplyr::contains("samp_"), ~ . * population),
                type = paste0("Number circumcised (",
                              stringr::str_remove(type, " coverage"),
                              ")")
            )
    })
    # also calculate unmet need
    tmp[[length(tmp) + 1]] <- results %>%
        dplyr::filter(type == "MC coverage") %>%
        dplyr::mutate(
            dplyr::across(dplyr::contains("samp_"), ~ population * (1 - .)),
            type = "Unmet need"
        )

    # Append together
    results_n <- as.data.frame(data.table::rbindlist(tmp, use.names = T))
}