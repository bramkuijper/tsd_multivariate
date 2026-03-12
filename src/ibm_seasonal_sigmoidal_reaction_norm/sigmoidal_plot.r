library("tidyverse")
library("purrr")

sigmoidal <- function(
       temp,
       tmax,
       tb,
       ta,
       ...)
{
    return(tmax / (1 + exp(-tb * (temp - ta))))
}


par.grid <- expand_grid(
    temp = seq(-5,5,.1),
    tmax = c(0.1,0.5,0.9),
    tb = seq(-1,1,1),
    t = seq(-1,1,1)
) %>% mutate(
    prob_breed = purrr::pmap_dbl(., \(temp, tmax, tb, t) {
        sigmoidal(temp, tmax, tb,t)
    }),
    tb_f = as_factor(tb)
)

ggplot(data = par.grid,
       mapping = aes(x = temp, y = prob_breed)) +
    geom_line(mapping = aes(colour = tb_f)) +
    facet_grid(tmax ~ t, labeller="label_both") +
    theme_classic()

ggsave(filename = "sigmoidal_plot.pdf")

temperature_function <- function(
        time_point, 
        season_duration) {
    sin(time_point * 2 * pi / season_duration)
}

# now make the sinusoidal environment
max_t <- 50                  
t_start_per_season <- 3
temp_data <- tibble(
    time = seq(0,200,1)
) %>% mutate(
    temp = temperature_function(time_point = time, 
                                season_duration = max_t
                                ),
)

t_xintercept <- seq(t_start_per_season, max(temp_data$time), by = max_t)

ggplot(data = temp_data,
       mapping = aes(x = time, y = temp)) +
    geom_line(colour = "purple", size=1.5) +
    geom_vline(xintercept = t_xintercept) +
    theme_classic()

