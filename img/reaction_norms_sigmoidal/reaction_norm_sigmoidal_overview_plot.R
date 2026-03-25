library("tidyverse")


the_data <- read.table(
        file="summary_all.csv", 
    sep=";", 
    header=T) %>% mutate(
        temp_plasticity = if_else(mu_tb > 0,"ON","OFF")
    ) %>% filter(
        time >= max(time) - 950
    )

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
    values_to = "fraction"
) %>% mutate(
    label = if_else(
        grepl(
            pattern="fem",juvenile_type),"Females","Males"),
    panel_label = if_else(
        mu_tb == 0, "No plasticity", "Plasticity"
    ),
    panel_label2 = if_else(
        mu_tb == 0, "A", "B"
    ),
    female_warmer = if_else(toptm <= toptf, "Females in warmer envts", "Males in warmer envts")
)

label_data <- data.frame(
    panel_label = rev(unique(the_data_l$panel_label)),
    identifier = LETTERS[1:length(unique(the_data_l$panel_label))]
)

subdat <- the_data_l %>% 
    filter(grepl(x = file, pattern="sim_seasonal_tsd_18_03_2026_152608_.*")) %>%
    filter(intercept_change == max(intercept_change))

ggplot(data = the_data_l, 
       mapping = aes(x = intercept_change / 0.707, y = fraction)) +
    geom_hline(yintercept = 1, colour = "lightgrey", linewidth = 0.75) +
    geom_point(mapping = aes(colour = label), alpha = 0.5) +
    scale_colour_brewer(palette = "Set1") +
    labs(x = "Change in temperature (% of 1 SD)",
         y = "Number of offspring after vs before perturbation") +
    facet_grid(female_warmer ~ panel_label) + 
    guides(colour = guide_legend(title = NULL)) +
    geom_text(size = 10, 
              data = label_data, 
              mapping = aes(x = max(the_data_l$intercept_change / 0.707), y = max(the_data_l$fraction), label = identifier)) +
    theme_classic(base_size = 16) +
    theme(strip.background = element_rect(colour = "transparent")
          , strip.text = element_text(size = 16, colour = "#000000"))  

ggsave(filename="main_ms_plot.pdf", width=12, height=8)