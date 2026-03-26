library("tidyverse")
library("patchwork")

# single simulations in which you see adaptation to change
# in a non-plastic vs plastic context

find.params <- function(filename) {

    f <- readLines(filename)

    seq.rev <- rev(seq(1,length(f),1))

    for (line_i in seq.rev)
    {
        if (length(grep("^\\d",f[[line_i]])) > 0)
        {
            return(line_i)
        }
    }
}

summary_data <- read.table("summary_sample_simulations.csv"
                           ,sep=";"
                           ,header=T)

subset_no_plasticity <- summary_data %>% filter(mu_tb == 0) %>% select(file)


file_normal_data <- basename(subset_no_plasticity[1,1])
distribution_file <- paste0("distribution_", file_normal_data)

param_line <- find.params(filename = file_normal_data)

data_no_plast <- read_delim(file=file_normal_data
        ,delim=";"
        ,n_max=param_line-1
        ,col_names=T)

# get the parameters
data_no_plast_params <- read_delim(file=file_normal_data
        ,delim=";"
        ,skip=param_line
        ,col_names=c("name","value")
        )

# put parameters into a better format
params <- data_no_plast_params$value
names(params) <- data_no_plast_params$name

# first no plasticity
data_no_plast_individuals <- read.table(file=distribution_file,
                                       sep=";",
                                       header=T)

# plot the temperature curve
sinusoidal <- function(time, season_time, time_change, intercept_change) {
    val <- sin(2 * pi * time / season_time) 
    val <- if_else(time > time_change, val + intercept_change, val)
    return(val)
}

max_time <- 100
max_t_season <- max_time / 2
# make the temperature time series
time_series_data <- tibble(
    time = seq(0,max_time),
) %>% mutate(
    temp = sinusoidal(time = time
                      ,season_time = max_t_season 
                      ,time_change = max_t_season
                      ,intercept_change = as.numeric(params["intercept_change"]))
)


# now draw reaction norms of a whole bunch of individuals
breeding_probability <- function(t, temp, threshold, tmax, ta, tb) {
    val <- tmax / (1.0 + exp(-tb*(temp - ta)))

    return(if_else(t %% max_t_season > threshold, 
                   (1.0 - val)^((t %% max_t_season) - threshold) * val, 0))
}

list_reaction_norms <- list()

max_individuals <- 1000

rows_latest_individuals <- 
    data_no_plast_individuals %>% filter(time == max(time))

# let's just plot the first 100
for (individual_idx in 1:max_individuals) {
    
    row <- rows_latest_individuals[individual_idx,]
    
    p_breed <- breeding_probability(t = time_series_data$time %% max_time, 
                         temp = time_series_data$temp, 
                         threshold = row$time_threshold, 
                         tmax = row$tmax, 
                         ta = row$temp_a, 
                         tb = row$temp_b)
    
    time_series_singular_individual <- data.frame(
        t = time_series_data$time,
        p_breed = p_breed,
        id = individual_idx
    )
    
    list_reaction_norms[[individual_idx]] <- time_series_singular_individual
}

all_time_series <- bind_rows(list_reaction_norms)


ggplot(data = time_series_data,
       mapping = aes(x = time, y = temp)) +
    geom_line(colour="darkblue") +
    geom_line(data = all_time_series,
              mapping = aes(x = t, y = p_breed, group = id), 
              linewidth = 0.25, alpha = 0.1) +
    labs(x = "Time", 
         y = "Breeding probability") + 
    theme_classic(base_size = 16)





