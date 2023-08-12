# fig to make the perturbation plot for dave's MS
library("tidyverse")
library("cowplot")

# obtain the final line of the actual simulation data of each file
get_last_line_number <- function(file_name)
{
    # get all the lines
    all_lines <- readLines(con=file_name)
    
    # obtain number of lines
    nlines <- length(all_lines)
    
    # search for last line from the end
    # (coz faster)
    for (line_idx in seq(nlines,1,-1))
    {
        # find the line that starts with a number using grep()
        # if grep returns something, anything, its length > 0
        # hence the return the index of that line and we are done
        if (length(grep(pattern="^\\d",all_lines[[line_idx]]))>0)
        {
            return(line_idx)
        }
    }
    
    # no line means no data hence quit
    stop("no line found")
} # end get_last_line_number()

generation.limits <- c(19000,23000)

prepare_data <- function(file_name)
{
    param_line <- get_last_line_number(
            file_name = file_name)

    data <- read_delim(file = file_name
            ,delim=";"
            ,n_max=param_line
            ) %>% filter(generation >= filter.vals[1] & generation <= filter.vals[2])

    return(data)
} # end prepare_data


get_data_per_time <- function(file.name)
{
    file.name <- as.character(file.name)
    last_line <- get_last_line_number(file.name)

    data.f <- read_delim(file=file.name
            ,delim=";"
            ,n_max=last_line - 1) |> filter(generation >= generation.limits[[1]] & generation <= generation.limits[[2]])

    return(data.f)
} # end function get_data_per_time()

get_sims <- function(data.subset)
{
    df_total <- NULL

    # loop through rows of dataframe and retrieve
    # I have looked at purrr::pmap, but after seeing
    # the utter shambolic documentation with all the poorly
    # explained .x and whatever
    # bullshit, 'no thank yo'
    for (row_idx in 1:nrow(data.subset))
    {
        sub_df <- get_data_per_time(file.name=data.subset[row_idx,"file"])

        sub_df$mu_df <- data.subset[row_idx,"mu_df"]
        sub_df$mu_sr <- data.subset[row_idx,"mu_sr"]
        sub_df$mu_b <- data.subset[row_idx,"mu_b"]
        sub_df$init_df <- data.subset[row_idx,"init_df"]
        sub_df$init_dm <- data.subset[row_idx,"init_dm"]
        sub_df$unique_id <- row_idx
        df_total <- bind_rows(df_total, sub_df)
    }

    return(df_total)
} # end get_sims


#### read in summary file
the.summ <- read_delim(
        file="summary_perturb_switch.csv"
        ,delim=";") %>% filter(s_pert2 == 0.1 & s_pert1 == 0.9 & init_df == 1.0 & init_dm == 1.0)

# now get all the data from the individual simulations
dat_all <- get_sims(the.summ) %>% mutate(
        unique_id_f = as_factor(unique_id)
        ,df = d1
        ,dm = d2) %>% arrange(generation)

mean.w <- mean(dat_all$wbar)
sd.w <- sd(dat_all$wbar)

# then calculate standardized fitness
dat_all <- mutate(dat_all
        ,wbar_std=(wbar - mean.w)/sd.w
        )

dat_tsd_only <- filter(dat_all, mu_df == 0 & mu_sr > 0 & mu_b == 0)

############### the actual plots ###############
w_range <- c(-5,1.5)
sr_range <- c(0,1)

tsd_only_w <- ggplot(data = dat_tsd_only,
        mapping = aes(x=generation
                ,y=wbar_std
                ,width=unique_id_f)) +
    geom_line(colour="#007bb3") +
    labs(x="Generation",y="Standardized growth rate",title="TSD") +
    theme_cowplot(12) +
    theme(plot.title=element_text(face="plain",hjust=0.5)) +
    ylim(w_range[1],w_range[2])

# sex ratio data
dat_tsd_only_long <- pivot_longer(dat_tsd_only
        ,cols=c("sr1","sr2","b","df","dm")
        ,names_to="trait"
        ,values_to="trait_value")

tsd_only_sr <- ggplot(data=dat_tsd_only_long,
        mapping = aes(x=generation
                ,y=trait_value
                ,colour=trait
                ,width=unique_id_f
                )) +
        geom_line() +
        labs(x="Generation",y="Sex ratio") +
        theme_cowplot(12) +
        ylim(sr_range[1], sr_range[2])

############# burrowing #############

dat_tsd_burrowing <- filter(dat_all,mu_df == 0 & mu_sr > 0 & mu_b > 0)

tsd_burrowing_w <- ggplot(data = dat_tsd_burrowing,
        mapping = aes(x=generation
                ,y=wbar_std
                ,width=unique_id_f)) +
    geom_line(colour="#007bb3") +
    labs(x="Generation",y=NULL, title="TSD + burrowing") +
    theme_cowplot(12) +
    theme(plot.title=element_text(face="plain",hjust=0.5)) +
    ylim(w_range[1],w_range[2])

# sex ratio data
dat_tsd_burrowing_long <- pivot_longer(dat_tsd_burrowing
        ,cols=c("sr1","sr2","b","df","dm")
        ,names_to="trait"
        ,values_to="trait_value")

tsd_burrowing_sr <- ggplot(data=dat_tsd_burrowing_long,
        mapping = aes(x=generation
                ,y=trait_value
                ,width=unique_id_f)) +
        geom_line(mapping=aes(colour=trait)) +
        labs(x="Generation",y=NULL) +
        theme_cowplot(12) +
        ylim(sr_range[1], sr_range[2])


############# dispersal #############

dat_tsd_dispersal <- filter(dat_all,mu_df > 0 & mu_sr > 0 & mu_b == 0)

tsd_dispersal_w <- ggplot(data = dat_tsd_dispersal,
        mapping = aes(x=generation
                ,y=wbar_std
                ,width=unique_id_f)) +
    geom_line(colour="#007bb3") +
    labs(x="Generation",y=NULL, title="TSD + dispersal") +
    theme_cowplot(12) +
    theme(plot.title=element_text(face="plain",hjust=0.5)) +
    ylim(w_range[1],w_range[2])

# sex ratio data
dat_tsd_dispersal_long <- pivot_longer(dat_tsd_dispersal
        ,cols=c("sr1","sr2","b","df","dm")
        ,names_to="trait"
        ,values_to="trait_value")

tsd_dispersal_sr <- ggplot(data=dat_tsd_dispersal_long,
        mapping = aes(x=generation
                ,y=trait_value
                ,width=unique_id_f)) +
        geom_line(mapping=aes(colour=trait)) +
        labs(x="Generation",y=NULL) +
        theme_cowplot(12) +
        ylim(sr_range[1], sr_range[2])

plot_grid(tsd_only_w
        ,tsd_only_sr
        ,tsd_burrowing_w
        ,tsd_burrowing_sr
        ,tsd_dispersal_w
        ,tsd_dispersal_sr
        ,nrow=2
        ,ncol=3
        ,byrow=F
        ,labels=LETTERS[1:6])

ggsave(file="plot_perturbation_switch.pdf",width=10)
