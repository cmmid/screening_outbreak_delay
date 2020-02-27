library(tidyverse)
library(data.table)
library(magrittr)
library(grid)
library(scales)

#allow for nesting in ggplot2::facet_wrap
source("./ggplot_nested_facets.R")

# function required for providing the shape and scale parameter from a Gamma
moment_match <- function(mean, var){
  list(shape = mean^2/var,
       rate  = mean/var,
       scale = var/mean) 
}

# given quantiles of a Gamma, what are its parameters?
source("gamma.parms.from.quantiles.R")

# function to determine time to onset of symptoms or severe disease from a Gamma
time_to_event <- function(n, mean, var){
  if (var > 0){
    parms <- moment_match(mean, var)
    return(rgamma(n, shape = parms$alpha, rate = parms$beta))
  } else{
    return(rep(mean, n))
  }
}

# Simulate incubation and infectious periods, recovery time and flight departure
# and arrival times
generate_histories <- function(input){
  
  with(input,
       data.frame(
         incu  = time_to_event(n = sims, mean =  mu_inc, var =  sigma_inc),
         inf   = time_to_event(sims, mu_inf, sigma_inf),
         recov = time_to_event(sims, mu_recov, sigma_recov)) %>%
         mutate(flight.departure = runif(sims, min = 0, max =2*(incu + inf)),
                flight.arrival   = flight.departure + dur_flight) )
  
}

# given infection histories above, what proportion of travellers end up being 
# caught at each step in the screening process?
calc_outcomes <- function(infection_histories, parms){
  
  infection_histories %>%
    mutate(
      hospitalised_prior_to_departure = inf + incu < flight.departure) %>%
    filter(hospitalised_prior_to_departure == FALSE) %>%
    mutate(
      exit_screening_label  = runif(n(), 0, 1) < parms$sens.exit /100,
      entry_screening_label = runif(n(), 0, 1) < parms$sens.entry/100) %>%
    mutate(symp_at_exit   = incu < flight.departure,
           symp_at_entry  = incu < flight.arrival,
           
           
           # EXIT SCREENING
           found_at_exit  = symp_at_exit &   exit_screening_label,
           missed_at_exit = symp_at_exit &  !exit_screening_label,
           
           
           
           # SEVERE SYMPTOMS REQUIRE HOSPITALISATION
           sev_at_exit    = 0,              # no hospitalised can exit
           sev_from_lat   = 
             (!symp_at_exit) &              # latent on exit
             (incu + inf < flight.arrival), # progress before arrival
           sev_from_symp  = 
             missed_at_exit &               # symp but not caught at exit
             (incu + inf < flight.arrival), # progress before arrival
           sev_at_entry   = sev_from_lat | sev_from_symp,
           
           found_at_entry = 
             symp_at_entry &                # symptomatic when entering
             !sev_at_entry &
             entry_screening_label ,        # found by screening
           
           found_at_entry_only = 
             found_at_entry &               # symptomatic when entering
             (!symp_at_exit)                # but latent when exiting
    )
  
}

# given the screening results, above, what proportion are detected at exit, at
# entry, or severe at entry?
calc_probs <- function(outcomes,
                       parms){
  
  
  outcomes %>%
    summarise(
      prop_sev_at_entry = (1-parms$prop.asy/100)*mean(sev_at_entry),
      
      
      prop_symp_at_exit = (1-parms$prop.asy/100)*mean(found_at_exit),
      
      prop_symp_at_entry = (1-parms$prop.asy/100)*mean(
        (missed_at_exit & found_at_entry & !sev_at_entry) |
          (found_at_entry_only & !sev_at_entry))
    ) %>%
    mutate(prop_undetected = 1 - (prop_sev_at_entry +
                                    prop_symp_at_exit +
                                    prop_symp_at_entry)) %>%
    as.list %>%
    return
  
}

# list of pathogens that may be worth considering as sensitivity
pathogen <- list(
  `nCoV-2019` = 
    data.frame(
      # (Li et al. (2020) NEJM)
      mu_inc    =  5.2,
      sigma_inc =  4.1,
      mu_inf    =  9.1,
      sigma_inf = 14.7,
      prop.asy  = 17,
      mu_recov  = 14.7,
      sigma_recov = 56.1
    ),
  
  # `nCoV-2019 (Backer)` = 
  #   data.frame(
  #     # (Backer et al., 2020)
  #     mu_inc    =  5.7,
  #     sigma_inc =  round(2.6^2,1),
  #     # (Huang et al., 2020)
  #     mu_inf    = 8.0,
  #     sigma_inf = round(((13-5)/1.35)^2,1), # using Higgins (2008) Cochrane Handbook
  #     prop.asy  = 17
  #   ),
  
  `SARS-like (2002)` = 
    data.frame(
      mu_inc    =  6.4,
      sigma_inc = 16.7,
      mu_inf    =  3.8,
      sigma_inf =  6.0,
      prop.asy  =  0.0),
  `Flu A/H1N1-like (2009)` =
    data.frame(
      mu_inc    =  4.3,
      sigma_inc =  1.05,
      mu_inf    =  9.3,
      sigma_inf =  0.7,
      prop.asy  = 16.0 # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4586318/ 
    ),
  `MERS-like (2012)` = 
    data.frame(
      mu_inc    =  5.5,
      sigma_inc = 6.25, # https://www.sciencedirect.com/science/article/pii/S1473309913703049?via%3Dihub#sec1
      mu_inf    =  5.0,  # https://www.nejm.org/doi/10.1056/NEJMoa1306742
      sigma_inf =  7.5,
      prop.asy  = 21.0  # https://doi.org/10.1016/j.tmaid.2018.12.003 citing https://www.who.int/csr/disease/coronavirus_infections/risk-assessment-august-2018.pdf?ua=1
    ),
  # for the app at http://cmmid-lshtm.shinyapps.io/traveller_screening
  Custom     = 
    data.frame(
      mu_inc    =  5.0,
      sigma_inc =  5.0,
      mu_inf    =  5.0,
      sigma_inf =  5.0,
      prop.asy  = 10.0)) %>%
  dplyr::bind_rows(., .id = "name")


# pretty preting some CI labels in the app
make_ci_label <- function(x){
  x <- round(x)
  return(sprintf("%i (%i, %i)", x[1], x[2], x[3]))
}


# function to generate travellers to have screening applied
generate_travellers <- function(input, i){
  tibble(i = i) %>% 
    mutate(probs=future_pmap(.f=generate_histories,
                             list(   
                               dur.flight = input$dur.flight,
                               dur_flight = input$dur.flight/24,
                               mu_inc     = input$mu_inc,
                               sigma_inc  = input$sigma_inc,
                               mu_inf     = input$mu_inf,
                               sigma_inf  = input$sigma_inf,
                               sens.exit  = input$sens.exit,
                               sens.entry = input$sens.entry,
                               prop.asy   = input$prop.asy,
                               sims       = i)) %>%
             calc_outcomes(., input) %>%
             calc_probs(., input)) %>% 
    unnest_wider(probs)
}

# function to take travellers and work out their detection probabilities
generate_probabilities <- function(travellers){
  travellers %>%
    pivot_longer(cols = c(prop_symp_at_exit,
                          prop_symp_at_entry,
                          prop_sev_at_entry,
                          prop_undetected),
                 names_to = "screening",
                 values_to = "prob") %>% 
    group_by(screening) %>% 
    summarise(mean_prob = mean(prob*100),
              lb_prob=quantile(probs=0.025,x=prob*100),
              ub_prob=quantile(probs=0.975,x=prob*100)) %>%  
    pivot_longer(cols=c(mean_prob,
                        lb_prob,ub_prob)) %>% 
    pivot_wider(names_from = screening, values_from = value)
}


# back-calculate variance from the lower and upper bands of a symmetric
# confidence interval of the mean
get_var <- function(lower, upper, n, alpha = 0.05){
  n*((lower - upper)/(2*qt(1 - alpha/2, n-1)))^2
}

# pretty labels for when travellers are detected with varying degrees of severity
assign_labels <- function(x){
  
  x %>% 
    mutate(label = case_when(
      found_at_exit       == TRUE ~     "detected at exit screening",
      sev_at_entry        == TRUE ~     "detected as severe on flight",
      found_at_entry_only == TRUE ~     "detected at entry screening",
      TRUE ~ "not detected"
    ))
  
}


# https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1003277
T0 <- function(R0, k, c = 0.5){

  # 1-c the probability of outbreak
  # c the probability of control
  
  -log(c)/log(R0)*(0.334 + 0.689/k + 0.408/R0 - 
                     0.507/(k*R0) - 0.356/(R0^2) + 0.467/(k*R0^2))
}

# pretty label for plotting
contact_labeller <- function(x){
  paste("Traveller sensitisation:",
        scales::percent(as.numeric(x)))
}

# nicer axis labels
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

# turning on and off to calculate theta
screening_masks <- 
  list(`Exit only` = c("not detected",
                       "detected at entry screening",
                       "detected as severe on flight"),
       `Exit and entry` = c("not detected"),
       `No screening` = c("not detected",
                          "detected at exit screening",
                          "detected at entry screening",
                          "detected as severe on flight")) 

# what proportion of arrivals go undetected?
calculate_theta <- function(outcomes,
                            screening_masks){
  
  map_df(.x = screening_masks,
         ~filter(outcomes %>% assign_labels,
                 label %in% .x),
         .id = "Screening") %>%
    count(Screening) %>%
    mutate(N = nrow(outcomes),
           frac_undetected = n/N) %>%
    select(-c(n, N)) %>%
    mutate(Screening = fct_rev(fct_inorder(Screening)))
  
}

# pretty labels for plotting
k_labels  <-
  list(`COVID-19`  = 0.54,
       `SARS-like` = 0.16,
       `Flu-like`  = 2) %>%
  map_df(~data.frame(k = .x),
         .id     = "k_label") %>%
  mutate(k_label = paste0("k = ",k, " (",k_label,")")) %>%
  mutate(k_label = fct_inorder(k_label))

# make plot where facet has lines varying by lambda
make_results_plot_arrival <- function(delays, x = NULL){
  
  delays %<>% mutate(`Traveller screening` = scales::percent(contact_tracing))
  
  palette <- RColorBrewer::brewer.pal(
    n = 1 + length(unique(delays$lambda)),
    "Reds")[-1]
  
  palette_interventions <- c("No screening" = "grey60",
                             "Exit only" = "lightskyblue",
                             "Exit and entry" = "coral")
  
  palette_sens <- RColorBrewer::brewer.pal(
    n = 1 + length(unique(delays$`Traveller screening`)),
    "Purples")[-1]
  
  if (!is.null(x)){
    delays <- inner_join(delays, x)
  }
  
  delays %<>% make_nesters
  
  delays %>%
    filter(delta_t >= 0.5) %>%
   
    ggplot(
      .
    ) + aes(
      x     = delta_t, 
      y     = 1-delta_t_ecdf, 
      color = factor(lambda),
      group = interaction(factor(lambda),`Screening`)
    ) +
    geom_line(size = 1) +
    scale_x_log10("Delays (days)",
                  expand = expand_scale(),
                  breaks = c(1, 3, 10, 30, 100),
                  labels = function(x){ifelse(x < 1, "", x)},
                  limits = c(0.5, 130)) +
    scale_color_manual(values = palette,
                       name = expression(lambda == Infected~travellers~per~week
                       )) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, by=0.25),
                       labels = percent) +
    ylab("Percentage of delays at least this long") + 
    annotation_logticks(sides="b") +
    theme_bw() + theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      strip.background = element_rect(color=NA)
    ) +
    facet_nested(
      #ensure same order as in k_labels object
      # put COVID-19 at top
      k_label + Screening ~ nester_rho + rho_f, 
      labeller=labeller(
        nester_rho    = label_value,
        rho_f = label_parsed
      )
    ) 
  
}

# make plot where facet has lines varying by rho
make_results_plot_sensitise <- function(delays, x = NULL){
  
  delays %<>% mutate(`Traveller screening` = scales::percent(contact_tracing))
  
  palette <- RColorBrewer::brewer.pal(
    n = 1 + length(unique(delays$lambda)),
    "Reds")[-1]
  
  palette_sens <- RColorBrewer::brewer.pal(
    n = 1 + length(unique(delays$`Traveller screening`)),
    "Greens")[-1]
  
  if (!is.null(x)){
    delays <- inner_join(delays, x)
  }
  
  delays %<>% make_nesters
  
  delays %>%
    filter(delta_t >= 0.5) %>%
    ggplot(
      .
    ) + aes(
      x     = delta_t, 
      y     = 1-delta_t_ecdf, 
      color = Screening,
      group = interaction(factor(lambda),`Screening`, contact_tracing)
    ) +
    geom_line(size = 1,
              aes(color = `Traveller sensitisation`)) +
    scale_x_log10(name   = "Delays (days)",
                  expand = expand_scale(),
                  breaks = c(1, 3, 10, 30, 100),
                  labels = function(x){ifelse(x < 1, "", x)},
                  limits = c(0.5, 130)) +
    scale_color_manual(values = palette_sens,
                       name = "Traveller sensitisation"
    ) +
    scale_y_continuous(name = "Percentage of delays at least this long",
                       limits = c(0, 1),
                       seq(0, 1, by=0.25),
                       labels = percent) +
    annotation_logticks(sides="b") +
    theme_bw() + theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      
      strip.background = element_rect(color=NA)
    ) +
    facet_nested(nester_lambda + k_label + lambda_f ~  Screening,
                 labeller = 
                   labeller(lambda_f = label_parsed,
                            nester_lambda = label_parsed))
  
}


# make plot where facet has lines varying by theta
make_results_plot_screening <- function(delays, x = NULL){
  
  palette <- RColorBrewer::brewer.pal(
    n = 1 + length(unique(delays$lambda)),
    "Reds")[-1]
  
  palette_interventions <- c("No screening"   = "grey80",
                             "Exit only"      = "#41b6c4",
                             "Exit and entry" = "#225ea8")
  if (!is.null(x)){
    delays <- inner_join(delays, x)
  }
  
  delays %<>% make_nesters
  
  
  delays %>%
    filter(delta_t >= 0.5) %>%
    ggplot(
      .
    ) + aes(
      x     = delta_t,
      y     = 1-delta_t_ecdf,
      color = Screening,
      group = interaction(Screening, `Traveller sensitisation`)
    ) +
    geom_line(size = 1) +
    scale_x_log10(name   = "Delays (days)",
                  expand = expand_scale(),
                  breaks = c(1, 3, 10, 30, 100),
                  labels = function(x){ifelse(x < 1, "", x)},
                  limits = c(0.5, 130)) +
    scale_color_manual(values = palette_interventions#,
                       #name = "Infected travellers per week"
    ) +
    scale_y_continuous(name = "Percentage of delays at least this long",
                       limits = c(0, 1),
                       seq(0, 1, by=0.25),
                       labels = percent) +
    annotation_logticks(sides="b") +
    theme_bw() + theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      strip.background = element_rect(color=NA)
    ) +
    facet_nested(nester_lambda + k_label + lambda_f  ~ nester_rho+ rho_f,
                 labeller = 
                   labeller(lambda_f = label_parsed,
                            nester_lambda = label_parsed,
                            nester_rho    = label_value,
                            rho_f = label_parsed)) 

}

# make the smaller plot for the manuscript
make_results_plot_screening_small <- function(delays, x = NULL){
  
  palette <- RColorBrewer::brewer.pal(
    n = 1 + length(unique(delays$lambda)),
    "Reds")[-1]
  
  
  palette_interventions <- c("No screening"   = "grey80",
                             "Exit only"      = "#41b6c4",
                             "Exit and entry" = "#225ea8")
  if (!is.null(x)){
    delays <- inner_join(delays, x)
  }
  
  delays %<>% make_nesters
  
  
  delays %>%
    filter(delta_t >= 0.5) %>%

    ggplot(
      .
    ) + aes(
      x     = delta_t,
      y     = 1-delta_t_ecdf,
      color = Screening,
      group = interaction(Screening, `Traveller sensitisation`)
    ) +
    geom_line(size = 1,
              aes(linetype = `Traveller sensitisation`)) +
    scale_x_log10(name   = "Delays (days)",
                  expand = expand_scale(),
                  breaks = c(1, 3, 10, 30, 100),
                  labels = function(x){ifelse(x < 1, "", x)},
                  limits = c(0.5, 130)) +
    scale_color_manual(values = palette_interventions#,
                       #name = "Infected travellers per week"
    ) +
    scale_y_continuous(name = "Percentage of delays at least this long",
                       limits = c(0, 1),
                       seq(0, 1, by=0.25),
                       labels = percent) +
    annotation_logticks(sides="b") +
    theme_bw() + theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      strip.background = element_rect(color=NA)
    ) +
    facet_nested(lambda_f ~ nester_lambda  ,
                 labeller = 
                   labeller(lambda_f = label_parsed,
                            nester_lambda = label_parsed,
                            nester_rho    = label_value,
                            rho_f = label_parsed)) +
    scale_linetype_manual(name = "Percentage of travellers sensitised",
                          values =
                            seq(
                              from = 2*length(levels(delays$rho_f)) - 1,
                              to =  1,
                              by = -2))
}


# interpolate the CCDF so that we aren't trying to plot millions of 
make_interp <- function(x){
  distinct(x, delta_t, delta_t_ecdf) %>%
    filter(is.finite(delta_t)) %>%
    with(., 
         approx(x = log10(delta_t),
                y = delta_t_ecdf,
                xout = log10(delta_t_out))) %>%
    data.frame(delta_t = 10^.$x,
               delta_t_ecdf = .$y) %>%
    select(-x, -y)
}

# this makes the variables required for facet_nested in the plot
make_nesters <- function(x){
  
  lambdas <- sort(unique(x$lambda))
  rhos    <- sort(unique(x$contact_tracing))
  
  mutate(x,
         nester_lambda = "lambda == infected~travellers~per~week",
         nester_rho = "Percentage of travellers sensitised",
         lambda_f = factor(lambda, 
                           levels = lambdas,
                           labels = paste("lambda ==", 
                                          lambdas)),
         `Traveller sensitisation` = scales::percent(contact_tracing),
         rho_f = factor(
           contact_tracing,
           levels = rhos,
           labels = paste("rho ==",
                          100*rhos,
                          "*'%'")),
         Screening = factor(Screening,
                            levels = c("No screening",
                                       "Exit only",
                                       "Exit and entry"),
                            ordered = T))
}

# in case you want to make a slightly less cluttered, but less tidy, table
unfill_vec <- function(x) {
  # https://github.com/tidyverse/tidyr/issues/250
  same <- x == dplyr::lag(x)
  ifelse(!is.na(same) & same, NA, x)
}


# make summary tables showing the quantiles of the relevant variable
probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
make_quantiles <- function(x, var){
  x %>%
    filter(!(contact_tracing == 0 & Screening == "No screening")) %>%
    group_by(Screening, contact_tracing, lambda, k) %>%
    nest %>%
    mutate(Quantiles = map(data, ~quantile(.x[[var]],
                                           probs)),
           Quantiles = map(Quantiles, ~bind_rows(.))) %>%
    unnest(Quantiles)
}

# for writing a table that contains the relevant quantiles
# requires a table of quantiles, above
make_results_table <- function(x){
  ungroup(x) %>%
    arrange(k, lambda, contact_tracing, Screening) %>%
    mutate(`Sensitisation` = paste0(100*contact_tracing, "%")) %>%
    select(k, lambda, Sensitisation, Screening, contains("%")) %>%
    # mutate_at(.vars = vars(lambda, Sensitisation, k, k_label),
    #           .funs = unfill_vec) %>%
    rename_at(.vars = vars(contains("%")),
              .funs = function(x){
                paste0(100 - parse_number(x), "%")
              }) %>%
    mutate_at(.vars = vars(contains("%")),
              .funs = round, digits = 0) %>%
    mutate_at(.vars = vars(contains("%")),
              .funs = function(x){gsub(pattern = "Inf", 
                                       replacement = "\u221E",
                                       x = x)})
}
