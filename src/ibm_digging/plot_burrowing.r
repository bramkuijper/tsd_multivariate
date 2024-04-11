library("tidyverse")
library("patchwork")

data_with <- read_delim(
    file="sim_seasonal_tsd_03_04_2024_222415_7"
    ,delim=";"
    ,n_max = 25000-1)

data_with <- data_with %>% mutate(id="with")

data_without <- read_delim(
    file="sim_seasonal_tsd_03_04_2024_225417_7"
    ,delim=";"
    ,n_max = 25000-1)

data_without <- data_without %>% mutate(id="without")

all_data <- bind_rows(data_with, data_without) %>%
    mutate(Actual_time = time / 50,
           ) %>%
    filter(Actual_time > 50000 & Actual_time < 150000) 

p_sd <- ggplot(data = all_data,
       mapping = aes(x = Actual_time, y = a)) +
    geom_line(mapping = aes(colour = id)) +
    ylab("ESD Threshold") +
    scale_colour_manual(
        values=c("#06b0fb","#ec2212"),
        labels=c("Digging evolving", "Digging constant")) +
    theme_classic() +
    theme(legend.title = element_blank()
          ,axis.text.x = element_blank()
          ,axis.title.x = element_blank())  +
    ggtitle("A")

p_time <- ggplot(data = all_data,
                 mapping = aes(x = Actual_time, y = depth)
                 ) +
    geom_line(mapping = aes(colour = id)) +
    ylab("Digging depth") +
    xlab("Time step") +
    theme_classic() +
    scale_colour_manual(
        values=c("#06b0fb","#ec2212"),
        labels=c("Digging evolving", "Digging constant")) +
    guides(colour="none") +
    theme(legend.title = element_blank()
          ,axis.text.x = element_blank()
          ,axis.title.x = element_blank())  +
    ggtitle("B")

all_data_l <- all_data %>% pivot_longer(
    cols = c(surviving_female_juvs, surviving_male_juvs),
    names_to = "type",
    values_to = "fitness"
)

all_data_l <- all_data_l %>% mutate(
    fitness_type=fct_cross(type, id)
)

p_fitness <- ggplot(data = all_data_l,
                 mapping = aes(x = Actual_time, y = fitness)
                 ) +
    geom_line(mapping = aes(colour = fitness_type)) +
    ylab("Fitness") +
    xlab("Time step") +
    scale_colour_manual(
        values=c("#06b0fb","#ec2212","#963a71","#70963a"),
        labels=c(
            "Daughters, digging evolving", 
            "Sons, digging evolving",
            "Daughters, digging constant",
            "Sons, digging constant"
            )) +
    theme_classic()  +
    ggtitle("C")

(p_sd / p_time / p_fitness)

ggsave(filename="plot_burrowing.pdf",width=8,height=8)