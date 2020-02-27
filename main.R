source("utils.R")
library(xlsx)

# define the parameters required for the simulation
input <- data.frame(dur.flight = 12, # hours 
                    sens.exit  = 86,
                    sens.entry = 86) %>%
  bind_cols(filter(pathogen, name == "nCoV-2019")) %>%
  mutate(dur_flight = dur.flight/24, # days
         sims = 1e5) %>%
  mutate(prop.asy = 17)

# SIMULATE TRAVELLERS
set.seed(2019)
histories <- generate_histories(input) %>%
  filter(inf + incu >= flight.departure) 

outcomes  <- data.table(calc_outcomes(histories, input))

set.seed(2019)
theta     <- calculate_theta(outcomes, screening_masks)

# SIMULATE OUTBREAKS

n_R0 <- 1e5

# R0_parms <- moment_match(2.2, 0.53) # variance from riou 
R0_parms <- gamma.parms.from.quantiles(q = c(1.4, 3.9),
                                       p = c(0.025, 0.975))

set.seed(2019)

delay_parameter_grid_new <- 
  
  ## PARAMETERS
  data.frame(R0 =
               rgamma(n     = n_R0,
                      shape = R0_parms$shape,
                      rate  = R0_parms$rate),
             rowid = 1:n_R0) %>%
  crossing(k_labels) %>%
  arrange(rowid) %>%
  mutate(T0 = T0(R0 = R0, k = k, c = 0.5),
         T0 = ifelse(T0 <= 0, Inf, T0)) %>%
  
  ## INTERVENTIONS AND ARRIVAL RATES
  crossing(lambda           = 10^seq(0,2),
           theta,
           contact_tracing  = c(0, 0.3, 0.5, 0.7)
  ) %>%
  
  ## DELAYS
  mutate(Rq      = (1 - contact_tracing)*R0,
         q0      = runif(n()),
         t0      = qgamma(q0, shape=T0, rate=lambda/7),
         Tq      = T0(Rq, k, 0.5),
         Tq      = ifelse(Tq <= 0, Inf, Tq),
         tdash   = qgamma(q0, shape=Tq, rate=frac_undetected*lambda/7),
         delta_t = tdash - t0,
         delta_t = ifelse(is.infinite(t0), Inf, delta_t)
    ) 

# Calculate summary statistics for relevant variables

delay_summaries     <- make_quantiles(delay_parameter_grid_new, "delta_t")
threshold_summaries <- make_quantiles(delay_parameter_grid_new, "T0") %>%
  group_by(k) %>%
  summarise_at(.vars = vars(contains("%")), .funs = median)

main_results <- expand.grid(contact_tracing = c(0, 0.5),
                            k = 0.54,
                            Screening = c("No screening",
                                          "Exit only",
                                          "Exit and entry")) %>%
  mutate(Screening = fct_inorder(Screening)) %>%
  arrange(contact_tracing, k, Screening) 

ungroup(delay_summaries) %>%
  inner_join(main_results) %>%
  make_results_table %>%
  write.xlsx2(x = ., file = "results/main_results.xlsx",
              showNA = FALSE)

delay_summaries %>%
  make_results_table %>%
  write.xlsx2(x = ., file = "results/all_results.xlsx",
              showNA = FALSE)

# What proportion of delays are finite?

delay_parameter_grid_new %>%
  filter(!(Screening == "No screening" & contact_tracing == 0)) %>%
  group_by(k, lambda, Screening, contact_tracing) %>%
  summarise(mean = round(mean(is.infinite(delta_t)), digits = 4)) %>% 
  ungroup %>%
  distinct(Screening,
           contact_tracing, 
           mean) %>%
  mutate(contact_tracing = scales::percent(contact_tracing)) %>%
  spread(contact_tracing, mean)

# Visualisation requires we calculate the ECDF and then interpolate

nesting_vars <- c("k_label", "k", "lambda",
                  "Screening", "contact_tracing")

delta_t_out   <- 10^seq(-1, 3, length.out = 101)


delays_to_plot <- delay_parameter_grid_new %>% 
  nest(data =  -all_of(nesting_vars))

delays_to_plot_processed <- 
  delays_to_plot %>%
  # filter out when there's no intervention
  filter(!(contact_tracing == 0 & Screening == "No screening")) %>%
  # split and unnest so we don't lose grouping variables
  rowwise() %>%
  group_split() %>%
  map(~unnest(.x, data)) %>%
  # calculate ecdf object and return value of ECDF at all delta_t
  map(~mutate(.x, delta_t_ecdf = ecdf(delta_t)(delta_t))) %>%
  # we can't interpolate at infinity, but we don't need to
  map(~filter(.x, is.finite(delta_t))) %>%
  # interpolate to reduce number of points being plotted
  map(~nest(.x, data =  -all_of(nesting_vars))) %>%
  map(~mutate(.x, data_interp = map(data, make_interp))) %>%
  map_df(~unnest(.x, data_interp)) %>%
  # allocate to quantile ranges in case we want to fill under line
  mutate(prob = cut(delta_t_ecdf,
                    breaks = c(0,probs,1),
                    include.lowest = T))


# generate plot in main manuscript
make_results_plot_screening_small(
  delays = delays_to_plot_processed,
  x      = main_results) %>%
  ggsave(plot = . ,
         filename = "results/alt_screening_small.pdf", height = 6.6, width = 4.4)

# generate plots from SI
plot_file_roots <- list(lambda = "make_results_plot_arrival",
                        rho    = "make_results_plot_sensitise",
                        screen = "make_results_plot_screening")

plot_file_roots %>%
  map(~do.call(what = .x, args= list(delays = delays_to_plot_processed)) %>%
        ggsave(filename = paste0("results/", .x, ".pdf"),
               width  = 13.2,
               height = 17.6,
               units = "in"))
