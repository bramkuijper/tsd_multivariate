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

subset_no_plasticity <- summary_data %>% filter(mu_tb == 0 & toptf > toptm) %>% select(file)
subset_plasticity <- summary_data %>% filter(mu_tb > 0 & toptf > toptm) %>% select(file)

file_normal_data <- basename(subset_no_plasticity[1,1])
file_normal_data_plastic <- basename(subset_plasticity[1,1])
distribution_file <- paste0("distribution_", file_normal_data)
distribution_file_plastic <- paste0("distribution_", file_normal_data_plastic)

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
                      ,time_change = max_t_season - 12
                      ,intercept_change = as.numeric(params["intercept_change"])),
    temp_pre = temp,
    temp_post = temp
)

time_series_data$temp_pre[seq(max_t_season - 12,nrow(time_series_data),1)]  <- NA
time_series_data$temp_post[seq(1,max_t_season - 12,1)]  <- NA

clip <- function(vec, min, max)
{
    vec[vec < min] <- min
    vec[vec > max] <- max
    
    return(vec)
}

# now draw reaction norms of a whole bunch of individuals
breeding_probability <- function(t, temp, threshold, tmax, ta, tb) {
    prob_breed <- tmax / (1.0 + exp(-tb*(temp - ta)))
    
    prob_breed <- clip(prob_breed, 0, 1)
    val <- numeric(length = length(t))
    
    for (i in 1:length(val))
    {
        if (i %% max_t_season < threshold) {
            val[[i]] <- 0
        } else {
            val[[i]] <- prob_breed[[i]]
            
            if (i %% max_t_season > 2)
            {
                # now look backwards to get previous values of prob_breed
                for (j in seq(1,(i %% max_t_season) - 1,1))
                {
                    if (j %% max_t_season < threshold) {
                        next
                    } else {
                        val[[i]] <- val[[i]] * (1.0 - prob_breed[[j %% max_t_season]])
                    }
                }
            }
        }
    }
    

    return(val)    
}

# function to retrieve time series of phenology from individual data
generate_time_series_sampled_individuals <- function(dataset, max_individuals=1000)
{
    list_reaction_norms <- list()
    
    rows_latest_individuals <- 
        dataset %>% filter(time == max(time))
    
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
    
    return(all_time_series)
}    

# first no plasticity
data_no_plast_individuals <- read.table(file=distribution_file,
                                       sep=";",
                                       header=T)
# then plasticity
data_plast_individuals <- read.table(file=distribution_file_plastic,
                                       sep=";",
                                       header=T)

all_time_series_unconditional <- 
    generate_time_series_sampled_individuals(dataset = data_no_plast_individuals)

all_time_series_conditional <- 
    generate_time_series_sampled_individuals(dataset = data_plast_individuals)

# ok now create optima for these values 
draw_optima <- function(temp, opt, omega)
{
    return(exp(-(temp - opt)^2/omega))
}

time_series_data <- time_series_data %>% mutate(
    optima_f = draw_optima(temp = temp, 
                           opt = as.numeric(params["toptf"]), 
                           omega = as.numeric(params["omegaf"])),
    optima_m = draw_optima(temp = temp, 
                           opt = as.numeric(params["toptm"]), 
                           omega = as.numeric(params["omegam"])),
)

color_pheno <- "#0476D9FF"
color_temperature <- "#01260EFF"
color_male_opt <- "#618C03FF"
color_female_opt <- "#BFB304FF"

plot_unconditional <- ggplot(data = time_series_data,
       mapping = aes(x = time, y = temp)) +
    geom_line(colour=color_temperature, mapping = aes(y = temp_pre)) +
    geom_line(colour=color_temperature, mapping = aes(y = temp_post)) +
    geom_line(colour=color_male_opt, mapping = aes(x =time, y = optima_f)) +
    geom_line(colour=color_female_opt, mapping = aes(x =time, y = optima_m)) +
    geom_line(data = all_time_series_unconditional,
              mapping = aes(x = t, y = p_breed, group = id), 
              linewidth = 0.25, 
              alpha = 0.2,
              colour = color_pheno) +
    labs(x = "", 
         y = "Environment",
         title = "A. No plasticity in breeding time") + 
    theme_classic(base_size = 16) +
    xlim(0,80) +
    scale_y_continuous(
        breaks = seq(-1,1.5,by = 0.5),
        sec.axis = sec_axis(~.,  
                               breaks = seq(0,1,by=0.25),
                               name = "Breeding probability")) +
    theme(
        axis.title.y.right = element_text(color=color_pheno),
        axis.text.y.right = element_text(color=color_pheno)
    ) 

plot_conditional <- ggplot(data = time_series_data,
       mapping = aes(x = time, y = temp)) +
    geom_line(colour=color_temperature, mapping = aes(y = temp_pre)) +
    geom_line(colour=color_temperature, mapping = aes(y = temp_post)) +
    geom_line(colour=color_male_opt, mapping = aes(x =time, y = optima_f)) +
    geom_line(colour=color_female_opt, mapping = aes(x =time, y = optima_m)) +
    geom_line(data = all_time_series_conditional,
              mapping = aes(x = t, y = p_breed, group = id), 
              linewidth = 0.25, 
              alpha = 0.2,
              colour = color_pheno) +
    labs(x = "Time", 
         y = "Environment",
         title = "B. Plasticity in breeding time") + 
    xlim(0,80) +
    scale_y_continuous(
        breaks = seq(-1,1.5,by = 0.5),
        sec.axis = sec_axis(~.,  
                               breaks = seq(0,1,by=0.25),
                               name = "Breeding probability")) +
    theme_classic(base_size = 16) +
    theme(
        axis.title.y.right = element_text(color=color_pheno),
        axis.text.y.right = element_text(color=color_pheno)
    ) 


plot_unconditional / plot_conditional

ggsave(filename = "plot_reaction_norm_phenology.pdf")