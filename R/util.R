# Various utility functions for functionality used on multiple occasions
# throughout package. Generally un-exported.

#### %||% ####

`%||%` <- function(x, y) { # no lint
  if (is.null(x)) y else x
}

#### read_circ_data (exported) ####

#' @title Function to read in Circumcision Data
#' @description Function to read in circumcision data to fit model. Handles
#' csv with \code{\link[data.table]{fread}} (but outputs data as a
#' `data.frame`), and geographical data with code{\link[sf]{read_sf}} (for which
#' it also adds unique identifiers for each `area_level`).
#' @param path Path to data.
#' @param filters Optional named vector, whose values dictate the values
#' filtered for in the corresponding column names. Only supports filtering for
#' one value for each column. default: NULL
#' @param selected Optional columns to select, removing others, default = NULL
#' @param ... Further arguments passed to or from other methods.
#' @seealso
#'  \code{\link[data.table]{fread}}
#'  \code{\link[sf]{read_sf}}
#' @return relevant data set, filtered as desired.
#' @rdname read_circ_data
#' @export
#' @importFrom dplyr %>%
#' @importFrom rlang .data
read_circ_data <- function(path, filters = NULL, selected = NULL, ...) {
  
  # maybe add a warning for missing "circ" columns for surveys?? And add
  # in NAs in this situation (look at KEN for this)
  
  # read in data, depending on file type
  cond <- tools::file_ext(path) %chin% c("geojson", "shp", "shx")
  if (cond) {
    .data <- sf::read_sf(path, ...)
  } else {
    # selection prior to loading is allowed by fread
    .data <- as.data.frame(data.table::fread(path, select = c(selected), ...))
  }
  
  # if desired, recursively filter data with provided `filters` vector
  if (!is.null(filters)) {
    cols <- names(filters)
    vals <- as.vector(filters[seq_along(filters)])
    for (i in seq_along(filters)) {
      if (!cols[i] %chin% names(.data)) next
      # change col i to symbol (if present), evaluate corresponding filter
      .data <- dplyr::filter(.data, !!rlang::sym(cols[i]) == vals[i])
    }
  }
  
  # Select specific columns, if desired (and present) (no need to do for fread)
  if (!is.null(selected) && cond) {
    .data <- .data %>%
      dplyr::select(dplyr::all_of(selected[selected %chin% names(.data)]))
  }
  
  # for areas, add unique identifier within Admin code and merge to boundaries
  if (inherits(.data, "sf")) {
    .data <- .data %>%
      dplyr::group_by(.data$area_level) %>%
      dplyr::mutate(space = dplyr::row_number()) %>%
      dplyr::ungroup()
  }
  return(.data)
}

#### create_dirs_r (exported) ####

#' @title Recursively Create Missing Directories
#' @description Function to recursively create directories if any of the
#' directories in a provided path are missing. Similar to \code{mkdir -p} from
#' Bash.
#' @param dir_path Path to a file or directory which you want to generate.
#' @rdname create_dirs_r
#' @export
create_dirs_r <- function(dir_path) {
  
  # split by "/"
  dirs <- stringr::str_split(dir_path, "/")[[1]]
  
  # check if dir_path is for dir (required) or specific file
  # does our string end in "/"? Then likely a dir
  cond <- substr(dir_path, nchar(dir_path), nchar(dir_path))
  cond <- grepl("/", cond)
  # does last word contain "."? Likely a file
  cond <- cond & !grepl(".", dplyr::last(dirs))
  if (!cond) {
    dirs <- dirs[-length(dirs)] # remove file name
  }
  if (length(dirs) == 0) stop("No missing directories created")
  # loop through each directory level in target directory, create any missing
  for (i in seq_along(dirs)) {
    # "recursively" check directories
    if (i == 1) {
      spec_dir <- dirs[i]
    } else {
      spec_dir <- paste(dirs[1:i], collapse = "/")
    }
    # create if missing
    if (!dir.exists(spec_dir)) {
      dir.create(spec_dir)
    }
  }
}

#### add_area_id ####

#' @title Change \code{area_id} from one hierarchy level to another
#' @description Function to change \code{area_id} from one hierarchy level to
#' another.
#' @param df Dataframe with \code{area_id} column.
#' @param df_areas_wide \code{sf} \code{dataframe} with shapefiles and area
#' hierarchy.
#' @param par list with two entries:
#' \itemize{
#'  \item{\code{area_lev}}{Current area level of \code{df}.}
#'  \item{\code{area_lev_select}}{Desired area level for \code{df}.}
#' }
#' @param add_keep_cols Additional columns to keep when summarising,
#' Default: NULL
#' @return \code{df} with `area_id` changed to `area_lev_select`.
#' @importFrom dplyr %>%
#' @rdname add_area_id
#' @keywords internal
add_area_id <- function(df,
                        df_areas_wide,
                        par,
                        add_keep_cols = NULL) {
  
  # Get area_id's
  area_lev_current_id <- paste0("area_id", par$area_lev)
  # The level we want
  area_lev_select_id <- paste0("area_id", par$area_lev_select)
  area_lev_select_name <- paste0("area_name", par$area_lev_select)
  
  # fix for when we don't want to aggregate fully to country level
  # choose consecutively higher area levels until we find one with shapefiles
  i <- 1
  area_lev_select_name_orig <- area_lev_select_name
  while (!area_lev_select_id %in% names(df_areas_wide)) {
    area_lev_select_id <- names(df_areas_wide)[i]
    area_lev_select_name <- paste0("area_name", i)
    i <- i + 1
  }
  if (i != 1) {
    paste0(
      "!(", 
      area_lev_select_name_orig, 
      " %in% names(df_areas_wide)), ", 
      "using ",
      area_lev_select_name, 
      " (i.e. aggregating to area level ", 
      i - 1,
      ")"
    )
  }
  
  # only select columns in our dataframe ("model" may be missing,
  # and `age` and `age_group` are interchangable)
  select_cols <- c("year", "age", "age_group", "population", "type", "model")
  # additional columns to keep, if supplied
  if (!is.null(add_keep_cols)) {
    select_cols <- unique(c(select_cols, add_keep_cols))
  }
  # remove columns which interfere with select below
  select_cols <- select_cols[!select_cols %chin% c("area_id", "area_name")]
  select_cols <- select_cols[select_cols %chin% names(df)]
  
  df_area_id <- df %>%
    # join in area names for chosen area_id
    dplyr::select(-dplyr::contains("area_name")) %>%
    dplyr::left_join(df_areas_wide %>%
                       dplyr::select(
                         # current level
                         area_id = dplyr::all_of(area_lev_current_id),
                         # desired level and
                         dplyr::all_of(area_lev_select_id),
                         # corresponding name
                         area_name = dplyr::all_of(area_lev_select_name)
                       ) %>%
                       dplyr::distinct(),
                     by = c("area_id")
    ) %>%
    # Select the right columns (account for when we are at the lowest level)
    dplyr::select(
      "area_id" = ifelse(par$area_lev_select == par$area_lev,
                         "area_id",
                         area_lev_select_id
      ),
      "area_name",
      dplyr::all_of(select_cols),
      dplyr::all_of(par$sample_cols)
    )
  
  # add missing area_level col, if required
  if (!"area_level" %in% names(df_area_id)) {
    df_area_id$area_level <- par$area_lev_select
  }
  
  return(df_area_id)
}


#### combine_areas ####

#' @title Collect results for lower area hierarchies
#' @description Function to collect results for lower area hierarchies by
#' joining higher area hierarchies.
#' @param .data Results for highest area hierarchy, to be combined to
#' give results for lower/less granular area hierarchies.
#' @param areas_wide \code{data.frame} with shapefiles and area hierarchy.
#' @param area_lev Desired admin boundary level.
#' @param join Indicator to decide whether to join data for different
#' area hierarchies, or return them in list form.
#' @param add_keep_cols Additional columns to keep when summarising,
#' Default: NULL
#' @param ... Further arguments passed to \link[threemc]{add_area_id}.
#' @return \code{data.frame} or list (depending on the value of \code{join})
#' with results for all area levels less than or equal to \code{area_lev}.
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @rdname combine_areas
#' @keywords internal
combine_areas <- function(.data,
                          areas_wide,
                          area_lev,
                          join,
                          add_keep_cols = NULL,
                          ...) {
  
  # add area_level to original df, if required
  if (!"area_level" %in% names(.data)) .data$area_level <- area_lev
  
  # all area levels in the data (0 indexed)
  area_levs <- seq_len(area_lev) - 1
  if (length(area_levs) == 0) area_levs <- -1
  
  # if only modelling a single area, no need to disaggregate!
  if (length(unique(.data$area_id)) == 1) area_levs <- area_lev
  
  # columns to keep
  add_keep_cols <- c(add_keep_cols, names(.data)[grepl("samp_", names(.data))])
  if (length(add_keep_cols) == 0) add_keep_cols <- NULL
  
  # collect results for lower area hierarchies by joining higher area
  # hierarchies (do so until you reach "area 0")
  if (area_levs[1] > -1) {
    results_list <- lapply(area_levs, function(x) {
      add_area_id(
        df = (.data %>% dplyr::select(-.data$area_name)),
        df_areas_wide = areas_wide,
        par = list(
          "area_lev" = area_lev,
          "area_lev_select" = x
        ),
        add_keep_cols = add_keep_cols
      )
    })
    # also add highest area hierarchy data
    results_list[[length(results_list) + 1]] <- .data
  } else {
    results_list <- .data
  }
  
  # return list or dataframe?
  if (join == TRUE) {
    return(data.table::rbindlist(results_list, use.names = TRUE, ...))
  }
  
  return(results_list)
}

#### append_mc_name ####

#' @title Append "mc" to Appropriate Column names
#' @description Ensure names for MC columns in fit have the suffix "_mc"
#' @param .data Dataframe/tibble whose columns include "MC" calculations.
#' @return \code{.data}, with column names appended appropriately.
#' @rdname append_mc_name
#' @keywords internal
append_mc_name <- function(.data) {
  mmc_tmc <- paste(c("mmc", "tmc"), collapse = "|")
  locs <- !(grepl(paste(c(mmc_tmc, "mc"), collapse = "|"), names(.data)))
  names(.data)[locs] <- paste0(names(.data)[locs], "_mc")
  
  return(.data)
}

#### create_aggregate_structure ####

#' @title Create a List containing the Hierarchy of levels
#'
#' @description Create a list containing all the area dependencies and number
#' of for each area in the hierarchy
#'
#' @inheritParams prepare_survey_data
#' @returns A list of length 2 containing:
#' \itemize{
#'  \item{"sub_region_list"}{A list of the specific sub-regions contained
#'  within each space (i.e. for each area_id) (including itself)}
#'  \item{"n_links_df"}{A dataframe with 2 columns, area_id and ndep, detailing
#'  the number of sub-regions contained within each area_id (also including
#'  itself)}
#' }
#' @rdname create_aggregate_structure
#' @keywords internal
create_aggregate_structure <- function(areas,
                                       area_lev) {
  # drop geometry and filter to specified area level
  areas <- sf::st_drop_geometry(areas) %>%
    dplyr::filter(.data$area_level <= area_lev)
  
  # Long to wide hierarchy (Need for new aggregation matrices)
  areas_wide <- spread_areas(areas)
  
  # iterate over all area_ids at specific area_lev (i.e. spaces)
  max_space <- areas %>%
    dplyr::filter(.data$area_level == area_lev) %>%
    dplyr::summarise(max(.data$space)) %>%
    dplyr::pull()

  if (max_space != nrow(areas)) {
    message(paste0(
        "max(areas$space) != nrow(areas), ",
        "may produce error in create_aggregate_structure"
    ))
  }

  area_id_seq <- seq(1, max_space, 1)
  sub_region_list <- lapply(area_id_seq, function(i) {
    # Get areas lower in the hierarchy
    areas_wide %>%
      dplyr::filter(dplyr::if_any(dplyr::starts_with("space"), ~ . == i)) %>%
      dplyr::pull(paste0("space", area_lev))
  })
  
  n_sub_region_df <- areas %>%
    dplyr::distinct(.data$area_id) %>%
    dplyr::mutate(sp_dep = vapply(sub_region_list, length, numeric(1)))
  
  return(
    list(sub_region_list = sub_region_list, n_sub_region_df = n_sub_region_df)
  )
}


#### spread_areas ####

#' Spread area hierarchy to wide format
#'
#' @inheritParams prepare_survey_data
#' @param min_level integer specifying the minimum level wanted
#' @param max_level integer specifying the maximum level wanted
#' @param space whether to include "space" columns. Excluding these returns the
#' same object as \code{naomi::spread_areas}, Default: TRUE
#'
#' @export
#' @importFrom dplyr %>%
#' @importFrom rlang .data :=
#' @rdname spread_areas
#' @keywords internal
spread_areas <- function(areas,
                         min_level = min(areas$area_level),
                         max_level = max(areas$area_level),
                         space = TRUE) {
  if (inherits(areas, "sf")) {
    boundaries <- areas %>% dplyr::select(.data$area_id)
    areas <- sf::st_drop_geometry(areas)
  } else {
    boundaries <- NULL
  }
  stopifnot(min_level >= min(areas$area_level))
  stopifnot(max_level <= max(areas$area_level))
  areas_wide <- areas %>%
    dplyr::filter(.data$area_level == min_level) %>%
    dplyr::select(
      `:=`(!!paste0("area_id", min_level), .data$area_id),
      `:=`(!!paste0("area_name", min_level), .data$area_name),
      `:=`(!!paste0("space", min_level), .data$space)
    )
  
  # create safe sequence from min_level + 1 to max_level to loop over
  level_seq <- seq_len(max_level)
  if (length(level_seq) == 0) {
    level_seq <- 0
  } else {
    level_seq <- level_seq[level_seq >= (min_level + 1)]
  }
  
  if (all(level_seq != 0)) {
    for (level in level_seq) {
      areas_wide <- areas_wide %>%
        dplyr::left_join(
          areas %>%
            dplyr::filter(.data$area_level == level) %>%
            dplyr::select(
              `:=`(!!paste0("area_id", level), .data$area_id),
              `:=`(!!paste0("area_name", level), .data$area_name),
              .data$parent_area_id,
              `:=`(!!paste0("space", level), .data$space),
            ),
          by = stats::setNames(
            c("parent_area_id"), c(paste0("area_id", level - 1L))
          )
        )
    }
  }
  
  areas_wide$area_id <- areas_wide[[paste0("area_id", max_level)]]
  if (!is.null(boundaries)) {
    areas_wide <- sf::st_as_sf(
      dplyr::left_join(areas_wide, boundaries, by = "area_id")
    )
  }
  
  # remove "space" columns returns same object as naomi::spread_areas
  if (space == FALSE) {
    areas_wide <- areas_wide %>%
      dplyr::select(-dplyr::contains("space"))
  }
  
  return(areas_wide)
}

#### match_age_group_to_ages ####

#' Create data frame of all ages within provided age group.
#'
#' @param age_group Age group, either "x-x" for a fixed upper age, or "x+", for
#' an age group with an upper age of `max_age`.
#' @param max_age Maximum age for age groups with no upper limit, Default: 60
#'
#' @rdname match_age_group_to_ages
#' @keywords internal
match_age_group_to_ages <- function(age_group, max_age = 60) {
  # if age_group ~ "x-x", expand age group from lower to upper age
  if (grepl("-", age_group)) {
    age_bounds <- as.integer(strsplit(age_group, "-")[[1]])
    ages <- age_bounds[1]:dplyr::last(age_bounds)
  } else {
    # if age group ~ "x+", take ages from x to max_age
    lower_age <- as.integer(gsub("+", "", age_group, fixed = TRUE))
    ages <- lower_age:max_age
  }
  # return data frame of single ages within provided age_group
  return(data.frame("age_group" = age_group, "age" = ages))
}

#' @title Change age group convention to match aggregation results
#' @description Change age group convention from "Y000_004" to "0-4", for
#' example.
#' @param .data `data.frame` with `age_group` column whose convention does not
#' match aggregations from `threemc`.
#' @rdname change_agegroup_convention
#' @keywords internal
change_agegroup_convention <- function(.data) {
  lower <- as.numeric(substr(.data$age_group, 3, 4))
  if (all(!is.na(as.numeric(lower)))) {
    upper <- as.numeric(substr(.data$age_group, 7, 8))
    .data$age_group <- paste(lower, upper, sep = "-")
  }
  return(.data)
}

#### survey_points_dmppt2_convert_convention ####

#' Create data frame of all ages within provided age group.
#' Convert survey coverage points & dmppt2 data to match convention of
#' aggregated results.
#'
#' @param .data Data frame with either survey calculated coverage, with
#' associated error bounds, or DMPPT2 coverage estimates calculated from VMMC
#' programme data.
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @importFrom data.table %chin%
#' @rdname survey_points_dmppt2_convert_convention
#' @export
survey_points_dmppt2_convert_convention <- function(.data) {
  
  # change column naming convention
  if ("survey_mid_calendar_quarter" %chin% names(.data)) {
    .data <- .data %>%
      dplyr::rename(
        year  = .data$survey_mid_calendar_quarter,
        type  = .data$indicator,
        mean  = .data$estimate,
        sd    = .data$std_error,
        lower = .data$ci_lower,
        upper = .data$ci_upper
      )
  } else if ("dmppt2_circumcision_coverage" %chin% names(.data)) {
    .data <- .data %>%
      dplyr::rename(mean = .data$dmppt2_circumcision_coverage)
  }
  # Change `type` column values to match that of threemc aggregations
  .data <- .data %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::matches("type"), ~ dplyr::case_when(
          . == "circumcised" ~ "MC coverage",
          . == "circ_medical" ~ "MMC coverage",
          TRUE ~ "TMC coverage"
        )
      )
    )
  
  # convert age group to threemc convention (i.e. Y000_004 -> 0-4)
  if (grepl("Y", .data$age_group[1])) {
    .data <- change_agegroup_convention(.data)
  }
  return(.data)
}

#### data.table internal functions ####

#' @title Convert vector to string
#' @description convert a vector like c(1, 4, 3, 2) into a string like
#' `[1, 4, 3, 2]` (common aggregation method for error messages). See also
#' `data.table:::brackify()`.
#' @rdname brackify
#' @keywords internal
brackify <- function(x, quote = FALSE) {
  cutoff <- 10L
  if (quote && is.character(x)) {
    x <- paste0("'", utils::head(x, cutoff + 1L), "'")
  }
  if (length(x) > cutoff) x <- c(x[seq_len(cutoff)], "...")
  sprintf("[%s]", toString(x))
}

#' @title Pattern matching for `data.table`
#' @description patterns returns the matching indices in the argument `cols`
#' corresponding to the regular expression patterns provided. The patterns must
#' be supported by `grep.` See also `data.table:::patterns()`.
#' @rdname patterns
#' @keywords internal
patterns <- function(..., cols = character(0L)) {
  l <- list(...)
  p <- unlist(l, use.names = any(nzchar(names(l))))
  if (!is.character(p)) stop("Input patterns must be of type character.")
  matched <- lapply(p, grep, cols)
  idx <- which(vapply(matched, length, FUN.VALUE = numeric(1)) == 0L)
  if (length(idx)) stop("Pattern(s) not found: [%s]", brackify(p[idx]))
  
  return(matched)
}

#### fill_downup_populations ####

#' @title Assume constant historical and future populations
#' @description Fills populations for historical/future years not present in 
#' populations dataset with earliest/latest known population for each unique 
#' area_id - age combination. 
#' @rdname fill_downup_populations
#' @keywords internal
fill_downup_populations <- function(
    populations, start_year, end_year, min_pop_year = NULL, max_pop_year = NULL
) {
  message(paste0(
    "Filling missing populations with earliest known population",
    " for each area_id and age"
  ))
  
  # calculate minimum and maximum year in populations if not provided
  if (is.null(min_pop_year)) min_pop_year <- min(populations$year)
  if (is.null(max_pop_year)) max_pop_year <- max(populations$year)
  
  # find years from provided start_year and min_pop_year
  missing_prev_years <- seq(start_year, min_pop_year - 1)
  missing_future_years <- seq(max_pop_year + 1, end_year)
  # add missing hostorical and future years to missing_years
  missing_years <- missing_prev_years
  if (missing_prev_years[2] < missing_prev_years[1]) {
    missing_years <- NULL
  }
  if (length(missing_future_years) > 1 && 
      (missing_future_years[2] < missing_future_years[1])) {
    missing_future_years <- NULL    
  }
  missing_years <- c(missing_years, missing_future_years)
  
  # create df matching populations columns for missing_years
  missing_rows <- tidyr::crossing(
    dplyr::select(populations, -c(.data$year, .data$population)),
    "year"       = missing_years,
    "population" = NA
  )
  # join to end of populations, order appropriately 
  populations <- dplyr::bind_rows(populations, missing_rows) %>%
    dplyr::arrange(.data$area_id, .data$age, .data$year) %>%
    # fill in with earliest known pop for each unique area_id & age combo
    dplyr::group_by(.data$area_id, .data$age) %>%
    tidyr::fill(.data$population, .direction = "downup") %>%
    dplyr::ungroup()
  
  return(populations)
}

#### naomi functions (not on CRAN) ####

#' @title Sample from TMB fit
#' @description See also `naomi::sample_tmb()`.
#' @rdname sample_tmb
#' @keywords internal
sample_tmb <- function(
    fit, nsample = 1000, rng_seed = NULL, random_only = TRUE, verbose = FALSE
  ) {
  
  set.seed(rng_seed)
  
  stopifnot(methods::is(fit, "naomi_fit"))
  stopifnot(nsample > 1)
  
  # internal function from TMB
  isNullPointer <- function(pointer) {
    .Call("isNullPointer", pointer) 
  }
  to_tape <- isNullPointer(fit$obj$env$ADFun$ptr)
  if (to_tape)
    fit$obj$retape(FALSE)
  
  if (!random_only) {
    if (verbose) print("Calculating joint precision")
    hess <- sdreport_joint_precision(fit$obj, fit$par.fixed)
    
    if (verbose) print("Inverting precision for joint covariance")
    cov <- solve(hess)
    
    if (verbose) print("Drawing sample")
    ## TODO: write version of rmvnorm that uses precision instead of covariance
    smp <- mvtnorm::rmvnorm(nsample, fit$par.full, cov)
    
  } else {
    r <- fit$obj$env$random
    par_f <- fit$par.full[-r]
    
    par_r <- fit$par.full[r]
    hess_r <- fit$obj$env$spHess(fit$par.full, random = TRUE)
    smp_r <- rmvnorm_sparseprec(nsample, par_r, hess_r)
    
    smp <- matrix(0, nsample, length(fit$par.full))
    smp[, r] <- smp_r
    smp[, -r] <- matrix(par_f, nsample, length(par_f), byrow = TRUE)
    colnames(smp)[r] <- colnames(smp_r)
    colnames(smp)[-r] <- names(par_f)
  }
  
  if (verbose) print("Simulating outputs")
  sim <- apply(smp, 1, fit$obj$report)
  
  r <- fit$obj$report()
  
  if (verbose) print("Returning sample")
  fit$sample <- Map(
    vapply, 
    list(sim), 
    "[[", 
    lapply(lengths(r), numeric), names(r)
  )
  is_vector <- vapply(fit$sample, inherits, logical(1), "numeric")
  
  fit$sample[is_vector] <- lapply(fit$sample[is_vector], matrix, nrow = 1)
  names(fit$sample) <- names(r)
  
  return(fit)
}

#' @title Get Joint Precision of TMB Fixed and Random
#' @description See also `naomi:::sdreport_joint_precision()`.
#' @rdname sdreport_joint_precision
#' @keywords internal
sdreport_joint_precision <- function(
    obj, 
    par.fixed               = NULL, 
    hessian.fixed           = NULL, 
    bias.correct            = FALSE, 
    bias.correct.control    = list(sd = FALSE, split = NULL, nsplit = NULL), 
    ignore.parm.uncertainty = FALSE, 
    skip.delta.method       = FALSE
  ) {
  
  if (is.null(obj$env$ADGrad) && (!is.null(obj$env$random))) {
    stop(paste0(
      "Cannot calculate sd's without type ADGrad available in object for ",
      "random effect models."
    ))
  }
  obj2 <- TMB::MakeADFun(
    obj$env$data, 
    obj$env$parameters, 
    type = "ADFun",
    ADreport = TRUE, 
    DLL = obj$env$DLL, 
    silent = obj$env$silent
  )
  r <- obj$env$random
  if (is.null(par.fixed)) {
    par <- obj$env$last.par.best
    if (!is.null(r)) {
      par.fixed <- par[-r]
    } else {
      par.fixed <- par
    }
    gradient.fixed <- obj$gr(par.fixed)
  } else {
    gradient.fixed <- obj$gr(par.fixed)
    par <- obj$env$last.par
  }
  if (length(par.fixed) == 0)
    ignore.parm.uncertainty <- TRUE
  if (ignore.parm.uncertainty) {
    hessian.fixed <- NULL
    pdHess <- TRUE
    Vtheta <- matrix(0, length(par.fixed), length(par.fixed))
  } else {
    if (is.null(hessian.fixed)) {
      hessian.fixed <- stats::optimHess(par.fixed, obj$fn, obj$gr)
    }
    pdHess <- !is.character(try(chol(hessian.fixed), silent = TRUE))
    Vtheta <- try(solve(hessian.fixed), silent = TRUE)
    if (methods::is(Vtheta, "try-error"))
      Vtheta <- hessian.fixed * NaN
  }
  if (!is.null(r)) {
    hessian.random <- obj$env$spHess(par, random = TRUE)
    L <- obj$env$L.created.by.newton
    if (!is.null(L)) {
      # non-exported function from TMB
      updateCholesky <- function(L, H, t = 0) {
        .Call("tmb_destructive_CHM_update", L, H, t)
      }
      updateCholesky(L, hessian.random)
      hessian.random@factors <- list(SPdCholesky = L)
    }
  }
  ADGradForward0Initialized <- FALSE
  ADGradForward0Initialize <- function() {
    obj$env$f(par, order = 0, type = "ADGrad")
    ADGradForward0Initialized <<- TRUE
  }
  if (!is.null(r)) {
    if (methods::is(L, "dCHMsuper")) {
      ## Non exported function from TMB
      solveSubset <- function(
        Q, 
        L = Matrix::Cholesky(Q, super = TRUE, perm = TRUE), diag = FALSE
      ) {
        stopifnot(methods::is(L, "dCHMsuper"))
        invQ <- .Call("tmb_invQ", L)
        iperm <- Matrix::invPerm(L@perm + 1L)
        if (diag) {
          invQ <- Matrix::diag(invQ)[iperm]
        } else {
          invQ <- invQ[iperm, iperm, drop = FALSE]
        }
        return(invQ)
      }
      diag.term1 <- solveSubset(L = L, diag = TRUE)
      if (ignore.parm.uncertainty) {
        diag.term2 <- 0
      } else {
        f <- obj$env$f
        w <- rep(0, length(par))
        if (!ADGradForward0Initialized) ADGradForward0Initialize()
        reverse.sweep <- function(i) {
          w[i] <- 1
          f(par, order = 1, type = "ADGrad", rangeweight = w,
            doforward = 0)[r]
        }
        nonr <- setdiff(seq_along(par), r)
        tmp <- sapply(nonr, reverse.sweep)
        if (!is.matrix(tmp))
          tmp <- matrix(tmp, ncol = length(nonr))
        A <- solve(hessian.random, tmp)
        diag.term2 <- rowSums((A %*% Vtheta) * A)
      }
      if (length(par.fixed) == 0) {
        jointPrecision <- hessian.random
      } else if (!ignore.parm.uncertainty) {
        G <- hessian.random %*% A
        G <- as.matrix(G)
        M1 <- methods::cbind2(hessian.random, G)
        M2 <- methods::cbind2(
          t(G), 
          as.matrix(t(A) %*% G) + hessian.fixed
        )
        M <- methods::rbind2(M1, M2)
        M <- Matrix::forceSymmetric(M, uplo = "L")
        dn <- c(names(par)[r], names(par[-r]))
        dimnames(M) <- list(dn, dn)
        p <- Matrix::invPerm(c(r, (seq_len(length(par)))[-r]))
        jointPrecision <- M[p, p]
      } else {
        message("ignore.parm.uncertainty ==> No joint precision available")
      }
    } else {
      message("Could not report sd's of full randomeffect vector.")
    }
  }
  return(jointPrecision)
}

#' @title Take multivariate normal sample from sparse precision matrix
#' @description See also `naomi:::rmvnorm_sparseprec()`.
#' @rdname rmvnorm_sparseprec
#' @keywords internal
rmvnorm_sparseprec <- function(
    n, 
    mean = rep(0, nrow(prec)), 
    prec = diag(length(mean))
  ) {
  
  z <- matrix(stats::rnorm(n * length(mean)), ncol = n)
  L_inv <- Matrix::Cholesky(prec)
  v <- mean + Matrix::solve(
    methods::as(L_inv, "pMatrix"), 
    Matrix::solve(Matrix::t(methods::as(L_inv, "Matrix")), z)
  )
  return(as.matrix(Matrix::t(v)))
}
