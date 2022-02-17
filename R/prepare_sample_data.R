# pull samples for various kinds of circumcision from TMB posterior fit
# and join with area data from "results"
prepare_sample_data <- function(N = 100,
                                populations,
                                no_prog_results = NULL,
                                prog_results = NULL,
                                no_prog_tmb_fit,
                                prog_tmb_fit,
                                type
) {

  if (is.null(no_prog_results) & is.null(prog_results)) {
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
      mmc <- "haz_mmc"
      tmc <- "haz_tmc"
      mc <- ifelse("haz" %in% names(tmp), "haz", "haz_mc")
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
      category = "coverage"
    } else category = type

    tmpx_1 <- tmpx_2 <- tmpx_3 <- tmpx_4 <- tmpx_5 <- tmpx_6 <- tmp
    if (type == "incidence") {
      tmpx_7 <- tmpx_8 <- tmpx_9 <- tmpx_10 <- tmpx_11 <- tmpx_12 <- tmp
    } else {
      tmpx_7 <- tmpx_8 <- tmpx_9 <- tmpx_10 <- tmpx_11 <- tmpx_12 <- NULL
    }

    #
    tmpx_1[, paste0("samp_", 1:N)] <- fit$sample[[mmc]][, 1:N]
    tmpx_1$type <- paste("MMC-nT", category)
    if (tmp$model[1] == "No program data") {

      tmpx_2[, paste0("samp_", 1:N)] <- 0
      tmpx_3[, paste0("samp_", 1:N)] <- fit$sample[[tmc]][, 1:N]
      tmpx_4[, paste0("samp_", 1:N)] <- fit$sample[[mmc]][, 1:N]
      tmpx_5[, paste0('samp_', 1:N)] <- fit$sample[[tmc]][, 1:N]
      # tmpx_6[, paste0("samp_", 1:N)] <- fit$sample[[tmc]][, 1:N] +
      #   fit$sample[[mmc]][, 1:N]
      tmpx_6[, paste0("samp_", 1:N)] <- fit$sample[[mc]][, 1:N]
    } else if (tmp$model[1] == "With program data") {

      if (type == "probability") {

        tmpx_2[, paste0("samp_", 1:N)] <- fit$sample$probs[, 1:N] *
          fit$sample[[tmc]][, 1:N]
        tmpx_3[, paste0("samp_", 1:N)] <- (1 - fit$sample$probs[, 1:N]) *
          fit$sample[[tmc]][ ,1:N]
        tmpx_4[, paste0("samp_", 1:N)] <- fit$sample[[mmc]][, 1:N] +
          fit$sample$probs[, 1:N] * fit$sample[[tmc]][, 1:N]
        tmpx_5[, paste0('samp_', 1:N)] <- fit$sample[[tmc]][, 1:N]
        # tmpx_6[, paste0("samp_", 1:N)] <- fit$sample[[tmc]][, 1:N] +
        #   fit$sample[[mmc]][, 1:N]
        tmpx_6[, paste0("samp_", 1:N)] <- fit$sample[[mc]][, 1:N]
      } else {

        tmpx_2[, paste0("samp_", 1:N)] <- fit$sample[[mmct]][,1:N]
        tmpx_3[, paste0("samp_", 1:N)] <- fit$sample[[tmc]][,1:N]
        tmpx_4[, paste0("samp_", 1:N)] <- fit$sample[[mmc]][,1:N] +
          fit$sample[[mmct]][,1:N]
        tmpx_5[, paste0("samp_", 1:N)] <- fit$sample[[tmc]][,1:N] +
          fit$sample[[mmct]][,1:N]
        # tmpx_6[, paste0("samp_", 1:N)] <- fit$sample[[mmc]][,1:N] +
        #   fit$sample[[tmc]][,1:N] + fit$sample[[mmct]][,1:N]
        tmpx_6[, paste0("samp_", 1:N)] <- fit$sample[[mc]][,1:N]
      }
    }
    tmpx_2$type <- paste("MMC-T", category)
    tmpx_3$type <- paste("TMC", category)
    tmpx_4$type <- paste("MMC", category)
    tmpx_5$type <- paste("TMIC", category)
    tmpx_6$type <- paste("MC", category)

    # Samples for the number of MCs performed (for incidence)
    if (type == "incidence") {
      tmpx_7 <- tmpx_1;  tmpx_7$type <- 'MMC-nTs performed'
      tmpx_8 <- tmpx_2;  tmpx_8$type <- 'MMC-Ts performed'
      tmpx_9 <- tmpx_3;  tmpx_9$type <- 'TMCs performed'
      tmpx_10 <- tmpx_4; tmpx_10$type <- 'MMCs performed'
      tmpx_11 <- tmpx_5; tmpx_11$type <- 'TMICs performed'
      tmpx_12 <- tmpx_6; tmpx_12$type <- 'MCs performed'
    }

    # Appending things together
    tmp <- as.list(mget(paste0("tmpx_", 1:12))) %>%
      bind_rows() %>%
      # only keep relevant columns
      select(area_id, area_name, year, age, type, model, contains("samp_")) %>%
      # join in region populations
      left_join(
        # only keep relevant columns in populations
        (populations %>%
           select(
             all_of(names(tmp)[names(tmp) %in% names(populations)]),
             population
           ))
      ) %>%
      relocate(population, .before = samp_1)

    # filter out na populations, with an appropriate message
    if (any(is.na(tmp$population)) == TRUE) {
        n1 <- nrow(tmp)
        tmp <- filter(tmp, !is.na(population))
        n2 <- nrow(tmp)
        if (n2 == 0) stop("No populations present in data")
        message(paste0("Missing population for ", n1 - n2, " records"))
    }


    return(tmp)
  }

  # Model with Probability of MC with no program data (only surveys)
  if (!is.null(no_prog_results)) {
    tmp1 <- no_prog_results %>%
      mutate(model = "No program data")
    tmp1 <- append_fun(tmp1, no_prog_tmb_fit, populations, type = type)
  } else {
    tmp1 <- NULL
  }

  # Model with Probability of MC with both programme and survey data
  if (!is.null(prog_results)) {
    tmp2 <- prog_results %>%
      mutate(model = "With program data")
    tmp2 <- append_fun(tmp2, prog_tmb_fit, populations, type = type)
  } else {
    tmp2 <- NULL
  }

  # Appending things together
  return(rbind(tmp1, tmp2))
}
