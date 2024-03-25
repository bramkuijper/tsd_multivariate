library("tidyverse")

temp <- seq(-1,1,0.01)

p_female <- function(a,b) {
    return(1.0 / (1.0 + exp(-(a + b*temp))))
}

p_female_x <- p_female(a=0,b=1.0)

df_data <- data.frame(temp=temp, pfem = p_female_x)

ggplot(data=df_data,
       mapping=aes(x=temp,y=pfem)) +
    geom_line() +
    theme_classic() +
    ylim(0,1)


# temperature 
temp <- seq(-1,1,0.01)

tf_opt <- -0.1
omega_f <- 0.1

# survival on females
psurv <- exp(-.5 * (tf_opt - temp)^2/omega_f)

df_surv <- 