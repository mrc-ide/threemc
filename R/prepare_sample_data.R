#' @title Pull N circumcision samples from TMB fit
#' @description Function to pull samples for various summaries inferred about
#' circumcision in a given area (prevalence, probability or incidence,
#' respectively).
#' @param N Number of samples to be generated, Default: 100
#' @param populations \code{data.frame} containing populations for each
#' region in tmb fits.
#' @param no_prog_results Quantiles of different summaries for different types
#' (Medical, traditional or total) circumcision, for survey data only,
#' Default: NULL
#' @param prog_results Quantiles of different summaries for different types
#' (Medical, traditional or total) circumcision, for survey data and VMMC
#' programme data, Default: NULL
#' @param no_prog_tmb_fit TMB model for Male Circumcision (MC), for survey data
#' only.
#' @param prog_tmb_fit TMB model for Male Circumcision (MC), for survey data and
#' VMMC programme data.
#' @param type Determines which aspect of MC in the regions in question you wish
#' to sample for. Can be one of "probability", "incidence" or "prevalence".
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
prepare_sample_data <- function(N = 100,
                                populations,
                                no_prog_results = NULL,
                                prog_results = NULL,
                                no_prog_tmb_fit,
                                prog_tmb_fit,
                                type) {
  if (is.null(no_prog_results) && is.null(prog_results)) {
    stop("cannot have prog_results == no_prog_results == NULL")
  }
  if (!type %in% c("probability", "incidence", "prevalence") |
    length(type) > 1) {
    stop("Please choose a valid type
         (one of 'probability', 'incidence', 'prevalence'")
  }

  if (nrow(populations) == 0) stop("No populations present in data")

  #########################################
  ### Loading rates from survival model ###
  #########################################

  # append samples from different circumcision models together (and pops)
  append_fun <- function(tmp, fit, populations, type) {

    # different objects to pull samples from fit, based on desired aggregation
    if (type == "probability") {
      mmc <- "haz_mmc" # medical circumcision
      tmc <- "haz_tmc" # traditional circumcision
      mc <- ifelse("haz" %in% names(tmp), "haz", "haz_mc") # all male circ
    } else if (type == "incidence") {
      mmc <- "inc_mmc"
      tmc <- "inc_tmc"
      mmct <- "inc_mmct"
      mc <- ifelse("inc" %in% names(tmp), "inc", "inc_mc")
    } else if (type == "prevalence") {
      mmc <- "cum_inc_mmc"
      tmc <- "cum_inc_tmc"
      mmct <- "cum_inc_mmct"
      mc <- ifelse("cum_inc" %in% names(tmp), "cum_inc", "cum_inc_mc")
    }

    # word to be pasted onto the end of circ type below
    if (type == "prevalence") {
      category <- "coverage"
    } else {
      category <- type
    }

    # initialise dataframes to store samples for different circ types
    tmpx_1 <- tmpx_2 <- tmpx_3 <- tmpx_4 <- tmpx_5 <- tmpx_6 <- tmp
    # type == "incidence" has 12 different "types"
    if (type == "incidence") {
      tmpx_7 <- tmpx_8 <- tmpx_9 <- tmpx_10 <- tmpx_11 <- tmpx_12 <- tmp
    } else {
      tmpx_7 <- tmpx_8 <- tmpx_9 <- tmpx_10 <- tmpx_11 <- tmpx_12 <- NULL
    }

    # Pull samples from data
    # Models with no VMMC data cannot distinguish between MMC-nT and MMC-T
    tmpx_1[, paste0("samp_", 1:N)] <- fit$sample[[mmc]][, 1:N]
    tmpx_1$type <- paste("MMC-nT", category)
    if (tmp$model[1] == "No program data") {
      tmpx_2[, paste0("samp_", 1:N)] <- 0
      tmpx_3[, paste0("samp_", 1:N)] <- fit$sample[[tmc]][, 1:N]
      tmpx_4[, paste0("samp_", 1:N)] <- fit$sample[[mmc]][, 1:N]
      tmpx_5[, paste0("samp_", 1:N)] <- fit$sample[[tmc]][, 1:N]
      tmpx_6[, paste0("samp_", 1:N)] <- fit$sample[[mc]][, 1:N]
    } else if (tmp$model[1] == "With program data") {
      if (type == "probability") {
        tmpx_2[, paste0("samp_", 1:N)] <- fit$sample$probs[, 1:N] *
          fit$sample[[tmc]][, 1:N]
        tmpx_3[, paste0("samp_", 1:N)] <- (1 - fit$sample$probs[, 1:N]) *
          fit$sample[[tmc]][, 1:N]
        tmpx_4[, paste0("samp_", 1:N)] <- fit$sample[[mmc]][, 1:N] +
          fit$sample$probs[, 1:N] * fit$sample[[tmc]][, 1:N]
        tmpx_5[, paste0("samp_", 1:N)] <- fit$sample[[tmc]][, 1:N]
        tmpx_6[, paste0("samp_", 1:N)] <- fit$sample[[mc]][, 1:N]
      } else {
        tmpx_2[, paste0("samp_", 1:N)] <- fit$sample[[mmct]][, 1:N]
        tmpx_3[, paste0("samp_", 1:N)] <- fit$sample[[tmc]][, 1:N]
        tmpx_4[, paste0("samp_", 1:N)] <- fit$sample[[mmc]][, 1:N] +
          fit$sample[[mmct]][, 1:N]
        tmpx_5[, paste0("samp_", 1:N)] <- fit$sample[[tmc]][, 1:N] +
          fit$sample[[mmct]][, 1:N]
        tmpx_6[, paste0("samp_", 1:N)] <- fit$sample[[mc]][, 1:N]
      }
    }
    # give appropriate labels to each df
    tmpx_2$type <- paste("MMC-T", category)
    tmpx_3$type <- paste("TMC", category)
    tmpx_4$type <- paste("MMC", category)
    tmpx_5$type <- paste("TMIC", category)
    tmpx_6$type <- paste("MC", category)

    # Samples for the number of MCs performed (for incidence)
    if (type == "incidence") {
      tmpx_7 <- tmpx_1
      tmpx_7$type <- "MMC-nTs performed"
      tmpx_8 <- tmpx_2
      tmpx_8$type <- "MMC-Ts performed"
      tmpx_9 <- tmpx_3
      tmpx_9$type <- "TMCs performed"
      tmpx_10 <- tmpx_4
      tmpx_10$type <- "MMCs performed"
      tmpx_11 <- tmpx_5
      tmpx_11$type <- "TMICs performed"
      tmpx_12 <- tmpx_6
      tmpx_12$type <- "MCs performed"
    }

    # Append together
    tmp <- as.list(mget(paste0("tmpx_", 1:12))) %>%
      dplyr::bind_rows() %>%
      # only keep relevant columns
      dplyr::select(
        .data$area_id, .data$area_name, 
        .data$year, .data$age, 
        .data$type, .data$model,
        dplyr::contains("samp_")
      ) %>%
      # join in region populations
      dplyr::left_join(
        # only keep relevant columns in populations
        (populations %>%
          dplyr::select(
            dplyr::all_of(names(tmp)[names(tmp) %in% names(populations)]),
            .data$population
          ))
      ) %>%
      dplyr::relocate(.data$population, .before = .data$samp_1)

    # filter out na populations, with an appropriate message
    if (any(is.na(tmp$population)) == TRUE) {
      n1 <- nrow(tmp)
      tmp <- tmp %>% dplyr::filter(!is.na(.data$population))
      n2 <- nrow(tmp)
      if (n2 == 0) stop("No populations present in data")
      message(paste0("Missing population for ", n1 - n2, " records"))
    }


    return(tmp)
  }

  # Model with Probability of MC with no program data (only surveys)
  if (!is.null(no_prog_results)) {
    tmp1 <- no_prog_results %>%
      dplyr::mutate(model = "No program data")
    tmp1 <- append_fun(tmp1, no_prog_tmb_fit, populations, type = type)
  } else {
    tmp1 <- NULL
  }

  # Model with Probability of MC with both programme and survey data
  if (!is.null(prog_results)) {
    tmp2 <- prog_results %>%
      dplyr::mutate(model = "With program data")
    tmp2 <- append_fun(tmp2, prog_tmb_fit, populations, type = type)
  } else {
    tmp2 <- NULL
  }

  # Append together
  return(rbind(tmp1, tmp2))
}
