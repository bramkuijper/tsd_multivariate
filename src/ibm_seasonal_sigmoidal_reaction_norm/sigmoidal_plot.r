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
    temp = seq(-1,1,1),
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
                  