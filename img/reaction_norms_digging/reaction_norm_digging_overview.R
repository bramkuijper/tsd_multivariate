library("tidyverse")

the_data <- read.table(file="summary_digging_plasticity.csv",sep=";",header=T)

# find out the different times pre and post perturbation
time_pre_post <- the_data %>% 
    dplyr::select(matches(match = "surviving_.*_juvs_")) %>%  # select column containing surviving_
    names() %>% # then get the column names of these columns
    gsub(pattern=".*?(\\d+)$",x=., replacement="\\1") %>%  # keep the numbers only
    as.numeric() %>% 
    unique() %>%
    sort()

the_data <- the_data %>% mutate(
   freq_fem_pre = !!rlang::sym(paste0("surviving_female_juvs_",time_pre_post[1])),
   freq_fem_post = !!rlang::sym(paste0("surviving_female_juvs_",time_pre_post[2])),
   freq_male_pre = !!rlang::sym(paste0("surviving_male_juvs_",time_pre_post[1])),
   freq_male_post = !!rlang::sym(paste0("surviving_male_juvs_",time_pre_post[2])),
   freq_fem_fraction = freq_fem_post / freq_fem_pre,
   freq_male_fraction = freq_male_post / freq_male_pre
)

the_data_l <- the_data %>% pivot_longer(
    cols = c(freq_fem_fraction, freq_male_fraction),
    names_to = "juvenile_type",
    values_to = "fraction") %>% mutate(
        female_warmer = if_else(toptm <= toptf, "Females in warmer envts", "Males in warmer envts"),
        is_female = if_else(grepl(x=juvenile_type,pattern = "fem"),"Female","Male")
    )

ggplot(data = the_data_l, 
       mapping = aes(x = intercept_change / 0.707, y = fraction)) +
    geom_hline(yintercept = 1, colour = "lightgrey", linewidth = 0.75) +
    geom_point(mapping = aes(colour = is_female), alpha = 0.5) +
    scale_colour_brewer(palette = "Set1") +
    facet_grid(female_warmer~mu_depth_slope) +
    theme_classic()
    
    
#+
#    scale_colour_brewer(palette = "Set1") +
#    labs(x = "Change in temperature (% of 1 SD)",
#         y = "Number of offspring after vs before perturbation") +
#    facet_grid(female_warmer ~ panel_label) + 
#    guides(colour = guide_legend(title = NULL)) +
#    geom_text(size = 10, 
#              data = label_data, 
#              mapping = aes(x = max(the_data_l$intercept_change / 0.707), y = max(the_data_l$fraction), label = identifier)) +
#    theme_classic(base_size = 16) +
#    theme(strip.background = element_rect(colour = "transparent")
#          , strip.text = element_text(size = 16, colour = "#000000"))  