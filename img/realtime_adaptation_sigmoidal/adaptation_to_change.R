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

params <- data_no_plast_params$value
names(params) <- data_no_plast_params$name

# first no plasticity
data_no_plast_individuals <- read.table(file=distribution_file,
                                       sep=";",
                                       header=T)

sinusoidal <- function(time, season_time, time_change, intercept_change) {
    val <- sin(2 * pi * time / season_time) 
    val <- if_else(time > time_change, val + intercept_change, val)
    return(val)
}

# make the temperature time series
probabilities_no_plast <- tibble(
    time = seq(0,100),
) %>% mutate(
    temp = sinusoidal(time = time
                      ,season_time = 50
                      ,time_change = 50
                      ,intercept_change = as.numeric(params["intercept_change"]))
)


# now draw reaction norms of a whole bunch of individuals
breeding_probability <- function(t, temp, threshold, tmax, ta, tb) {
    val <- tmax / (1.0 + exp(-tb*(temp - ta)))

    return(if_else(t >= threshold, val, 0))
}

list_reaction_norms <- list()

max_individuals <- 100

rows_latest_individuals <- 
    data_no_plast_individuals %>% filter(time == max(time))

# let's just plot the first 100
for (individual_idx in 1:max_individuals) {
    
    row <- rows_latest_individuals[individual_idx,]
    
    list_reaction_norms[[individual_idx ]]
    
}


ggplot(data = probabilities_no_plast,
       mapping = aes(x = time, y = temp)) +
    geom_line(colour="darkblue") +
    labs(x = "Time", 
         y = "Breeding probability") + 
    theme_classic(base_size = 16)





