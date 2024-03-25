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
tm_opt <- 0.0
omega_f <- 0.1
omega_m <- 0.1

# survival on females
psurv <- function(t_opt, omega, t)
{
    exp(-.5 * (t_opt - t)^2/omega)
}

df_surv <- data.frame(temp = temp, 
                      psurvf = psurv(t_opt = tf_opt,
                                     omega = omega_f,
                                     t = temp),
                      psurvm = psurv(t_opt = tm_opt,
                                     omega = omega_m,
                                     t = temp)
)

df_surv_l <- pivot_longer(data = df_surv,
                          cols=c(psurvf,psurvm),
                          names_to="Sex",
                          values_to="Survival"
                          )
                      
ggplot(data=df_surv_l,
       mapping=aes(x=temp, y = psurv)) +
    geom_line() +
    theme_classic()