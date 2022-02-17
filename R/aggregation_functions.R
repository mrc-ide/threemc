#  code common amongst `aggregationCode` sheets

# aggregate by area, year and type, for each individual age group modelled
# (weighted by population), and then convert to a percentage/probability
aggregate_sample_age_group <- function(
    results_list,
    add_groups,
    remove_groups) {

    if(inherits(results_list, "data.frame")) {
        stop("requires list from combine_areas (set argument join = FALSE)")
    }

    # standard age groups
    age_groups <- c('0-4',   '5-9',   '10-14', '15-19', '20-24', '25-29',
                    '30-34', '35-39', '40-44', '45-49', '50-54', '54-59',
                    '0+',    '10+',   '15+',   '15-24', '10-24',
                    '15-29', '10-29', '15-39', '10-39', '15-49', '10-49')
    # add or remove groups as desired
    if (!missing(add_groups)) age_groups <- c(age_groups, add_groups)
    if (!missing(remove_groups)) {
        age_groups <- age_groups[-(age_groups %in% remove_groups)]
    }

    # Multiplying by population to population weight
    results_list <- lapply(results_list, function(x) {
        x %>%
            mutate(across(contains("samp_"), ~ . * population))
    })
    # aggregate sample for each age group
    results <- lapply(seq_along(age_groups), function(i) {
        # If upper limit use this split
        if (grepl('-', age_groups[i]) == TRUE) {
            age1 <- as.numeric(strsplit(age_groups[i], '-')[[1]][1])
            age2 <- as.numeric(strsplit(age_groups[i], '-')[[1]][2])
        }
        # If no upper limit use this split
        if (grepl('\\+', age_groups[i]) == TRUE) {
            age1 <-  as.numeric(strsplit(age_groups[i], '\\+')[[1]][1])
            age2 <- Inf
        }
        results_list_loop <- lapply(results_list, function(x) {
            x <- x %>%
                # take results for age group i
                filter(age >= age1, age <= age2) %>%
                select(-age)
            # Getting summarising samples
            # group_by(area_id, area_name, year, model, type) %>%
            # summarise(across(everything(), sum), .groups = "drop") %>%

            x <- setDT(x)
            x <- x[,
                   lapply(.SD, sum, na.rm = T),
                   by = c("area_id", "area_name", "year", "model", "type"),
                   .SDcols = c("population", paste0("samp_", c(1 : N)))
            ]
            x <- x %>%
                # Adding age group
                mutate(age_group = age_groups[i])
        })
        # Printing index
        print(age_groups[i])
        # return ages
        return(results_list_loop)
        # Appending together
        # results <- rbind(results, rbindlist(results_list_loop))
    })
    # join together
    results <- as.data.frame(rbindlist(lapply(results, rbindlist)))

    # Multiplying by population to population weight
    # (don't do this for "N performed", if present)
    results <- results %>%
        mutate(across(contains("samp_"), ~ ifelse(grepl("performed", type),
                                                  .,
                                                  . / population))
        )


    return(results)
}

# aggregate survey points for each type and age group
aggregate_sample_survey <- function(
    survey_circumcision,
    join = TRUE,
    # types list
    types = list("Total" = c(unique(survey_circumcision$type)),
                 "Medical" = "MMC",
                 "Traditional" = "TMC"),
    age_groups = c("0-4",   "5-9",   "10-14", "15-19", "20-24", "25-29",
                   "30-34", "35-39", "40-44", "45-49", "50-54", "54-59",
                   "60-64", "65+",   "0+",    "10+",   "15+",   "15-24",
                   "15-29", "15-39", "15-49", "10-29", "10-39", "10-49",
                   "10-24")
) {

    # survey years
    survey_years <- unique(survey_circumcision$year)

    # loop through each type
    results_surv <- lapply(seq_along(types), function(i) {

        print(names(types)[i])

        # loop for each age group in the data (for each type)
        results_surv_type <- lapply(seq_along(age_groups), function(j) {

            print(age_groups[j])

            if (grepl("-", age_groups[j]) == TRUE) {
                age1 <- as.numeric(strsplit(age_groups[j], "-")[[1]][1])
                age2 <- as.numeric(strsplit(age_groups[j], "-")[[1]][2])
            }
            # If no upper limit use this split
            if (grepl("\\+", age_groups[j]) == TRUE) {
                age1 <-  as.numeric(strsplit(age_groups[j], "+")[[1]][1])
                age2 <- Inf
            }
            # Getting proportions
            tmp <- survey_circumcision %>%
                filter(age >= age1, age <= age2) %>%
                group_by(area_id, year) %>%
                summarise(
                    Y_ind =  length(circ_status) *
                        sum((circ_status == 1 & type %in% types[[i]]) *
                                indweight, na.rm = TRUE) /
                        sum(indweight, na.rm = TRUE),
                    Y_obs = sum(circ_status == 1 & type %in% types[[i]]),
                    N = length(circ_status),
                    p_ind = Y_ind / N,
                    p_obs = Y_obs / N,
                    .groups = "drop"
                ) %>%
                # adding age group
                mutate(age_group = age_groups[i])
            return(tmp)
        })
        # append together results for each age group
        results_surv_type <- as.data.frame(
            rbindlist(results_surv_type, use.names = T)
        )

        # Adding to skeleton dataset and adding regional information
        results_surv_type <- expand.grid(
            area_id = sort(unique(results_surv_type$area_id)),
            year = survey_years,
            type = names(types)[i],
            age_group = age_groups
        ) %>%
            left_join(results_surv_type) %>%
            # Adding region information
            left_join(
                (areas %>%
                     st_drop_geometry() %>%
                     select(contains("area"), -area_level_label)),
                by = "area_id"
            ) %>%
            left_join(
                (areas %>%
                     st_drop_geometry() %>%
                     select(
                         parent_area_id = area_id,
                         parent_area_name = area_name
                     )),
                by = "parent_area_id"
            )
        return(results_surv_type)
    })

    if (join == TRUE) results_surv <- as.data.frame(rbindlist(results_surv))

    return(results_surv)
}


# Getting summary statistics
posterior_summary_fun <- function(.data, N = 100) {

    # ensure numeric columns are after categorical
    .data <- .data %>%
        relocate(population | contains("samp_"), .after = everything())

    # pull locations of columns to "group by"
    id_cols <- seq_along(names(.data)[!names(.data) %like% "samp"])

    # use data.table as this can be quite slow for larger countries
    if(!inherits(.data, "data.table")) .data <- setDT(.data)

    # pivot to calculate row means and sds of samples for each stratification
    .data_long <- melt(.data, id.vars = id_cols,
                       measure.vars = c(paste0("samp_", 1:100)))
    .data <- .data_long[,
               '.'
               (mean = mean(value, na.rm = TRUE),
                   sd = sd(value, na.rm = TRUE)),
               keyby = c(names(.data)[id_cols])] # group by all categories]

    # calculate median and CI
    quantiles <- .data_long[, {
        quantiles = quantile(value, c(0.5, 0.025, 0.975), na.rm = TRUE, names = FALSE)
        ':='
        list(
            median = quantiles[1],
            lower = quantiles[2],
            upper = quantiles[3]
        )
    }, keyby = c(names(.data)[id_cols]), ]

    .data <- as.data.frame(merge(.data, quantiles))
}

# function to get change in prevalence/coverage from a given year
prevalence_change <- function(results, spec_year) {
    # pull samples from coverage in chosen year
    spec_year_results <- results %>%
        filter(year == spec_year) %>%
        select(-c(year, population)) %>%
        tidyr::pivot_longer(contains("samp_"), values_to = "prev_value")

    # join into results_change_2008 for corresponding categorical variables and subtract
    results_change_year <- results %>%
        tidyr::pivot_longer(contains("samp_")) %>%
        left_join(spec_year_results) %>%
        mutate(value = value - prev_value) %>%
        select(-prev_value) %>%
        tidyr::pivot_wider(., names_from = name, values_from = value) %>%
        mutate(type = paste0("Change in ", type, " from 2008"))
}

# calculate number of people circumcised (as well as unmet need)
n_circumcised <- function(results) {
    # Getting number of circumcised men
    tmp <- split(results, results$type)

    # get circumcised population by type
    tmp <- lapply(tmp, function(x) {
        x %>%
            mutate(
                across(contains("samp_"), ~ . * population),
                type = paste0("Number circumcised (",
                              stringr::str_remove(type, " coverage"),
                              ")")
            )
    })
    # also calculate unmet need
    tmp[[length(tmp) + 1]] <- results %>%
        filter(type == "MC coverage") %>%
        mutate(
            across(contains("samp_"), ~ population * (1 - .)),
            type = "Unmet need"
        )

    # Append together
    results_n <- as.data.frame(data.table::rbindlist(tmp, use.names = T))
}

# Merge regional information on the dataset (i.e. parent area info)
merge_area_info <- function(results, areas) {

    # Merging regional information on the dataset (i.e. parent area info)
    results <- results %>%
        # Adding region information
        left_join(
            (areas %>%
                 select(area_id:area_level)),
            by = c("area_id", "area_name")
        ) %>%
        relocate("area_level", .after = "area_name") %>%
        left_join(
            (areas %>%
                 select(parent_area_id = area_id,
                        parent_area_name = area_name)),
            by = "parent_area_id"
        ) %>%
        relocate(contains("area"))
}
