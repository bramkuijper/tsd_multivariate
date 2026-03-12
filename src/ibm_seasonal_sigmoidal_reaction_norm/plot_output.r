#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("tidyverse", warn.conflicts=F))
suppressPackageStartupMessages(library("jsonlite", warn.conflicts=F))
suppressPackageStartupMessages(library("patchwork", warn.conflicts=F))

# from a list of values like x1, x2, x3
# create a reasonable variable name, like x
make.var.name <- function(vars) {
    
    var1 <- vars[[1]]

    return(gsub(pattern="[_0-1]",replacement="",x=var1))
}

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
# use tidyverse's sym() to transform strings to symbols 
transform.sym  <- function(x) 
{
    if (!is.null(x))
    {
        sym(x)
    }
}

jsonstuff <- '[
    {"xvar" : "time",
        "yvar" : "a"
    },
    {
        "xvar" : "time",
        "yvar" : "b"
    },
    {
        "xvar" : "time",
        "yvar" : "temp_a"
    },
    {
        "xvar" : "time",
        "yvar" : "temp_b"
    },
    {
        "xvar" : "time",
        "yvar" : "tmax"
    },
    {
        "xvar" : "time",
        "yvar" : "time_threshold"
    },
    {
        "xvar" : "time",
        "yvar" : ["var_a","var_b","var_temp_a","var_temp_b","var_tmax","var_time_threshold"]
    },
    {
        "xvar" : "time",
        "yvar" : ["surviving_female_juvs","surviving_male_juvs"]
    },
    {
        "xvar" : "time",
        "yvar" : ["surviving_female_adults","surviving_male_adults"]
    },
    {
        "xvar" : "time",
        "yvar" : ["surviving_juv_sr","adult_sr","adult_surv_sr"]
    },
    {
        "xvar" : "time",
        "yvar" : ["available_males","available_females"]
    },
    {
        "xvar" : "time",
        "yvar" : ["already_attempted_males","already_attempted_females"]
    },
    {
        "xvar" : "time",
        "yvar" : ["nf","nm"]
    },
    {
        "xvar" : "time",
        "yvar" : "mean_environment"
    },
    {
        "xvar" : "time",
        "yvar" : "var_environment"
    }
]
'

distribution_file_prefix = "distribution_"

if (!exists("file.name"))
{
    
    # get command line arguments
    args = commandArgs(trailingOnly=TRUE)
    
    # give an error message if you do not provide it with a simulation file name
    if (length(args) < 1)
    {
        print("provide a simulation file name")
        stop()
    }
    
    file.name <- args[[1]]
}

param.line <- find.params(file.name)

data.tibble <- read_delim(file=file.name
        ,delim=";"
        ,n_max=param.line-1
        ,col_names=T)

# get the parameters
data.tibble.params <- read_delim(file=file.name
        ,delim=";"
        ,skip=param.line
        ,col_names=c("name","value")
        )

# transpose the tibble with the parameters
params <- data.tibble.params %>% pivot_wider(
        names_from = name
        ,values_from = value)

plot.structure <- fromJSON(jsonstuff, simplifyVector = F)

plot.structure.l <- length(plot.structure)

# list with all the plots
plot.list <- list()

plot.list.idx <- 1

data.tibble.middle <- data.tibble %>% mutate(
    time_diff = lead(time) - time,
    time_diff2 = abs(lag(time) - time)
) %>% filter(time_diff == 1 | time_diff2 == 1)

single.plot <- function(xvar, yvar, sub.data.tibble) {
    
    if (length(yvar) > 1)
    {
        yvar_name <- make.var.name(yvar)
        yvar_values <- paste0(yvar_name,"_values")

        sub.data <- pivot_longer(data=sub.data.tibble
                ,cols=yvar
                ,names_to=yvar_name
                ,values_to=yvar_values)

        # get rid of aes_string like this: 
        # https://stackoverflow.com/questions/74414272/how-to-replace-the-deprecated-ggplot2-function-aes-string-accepting-an-arbitrar/74414389#74414389 

        # aes arguments for the ggplot() call
        plot_args <- lapply(X=list(
                        x=xvar,
                        y=yvar_values),
                        FUN=transform.sym)

        # aes arguments for the geom_line() call
        line_args <- lapply(X=list(
                        colour=yvar_name),
                FUN=transform.sym)

        return(ggplot(data=sub.data
                ,mapping=aes(!!!plot_args)) + 
                    geom_line(mapping=aes(!!!line_args)))
    } else {
        
        # aes arguments for the ggplot() call
        plot_args <- lapply(X=list(
                        x=xvar,
                        y=yvar
                        ),
                FUN=transform.sym)

        return(ggplot(data=sub.data.tibble
                ,mapping=aes(!!!plot_args)) + geom_line())
    }
} # end single.plot

# first the 'normal' plots over the whole of
# evolutionary time
# then later on during the more 'dense' data sampling
# period around the change point
for (plot_struct_idx in 1:plot.structure.l)
{
    # get the (potential list of) y variable(s)
    # as this is a list and hence highly structured
    # hence, try to flatten it
    xvar <- unlist(plot.structure[[plot_struct_idx]]$xvar)
    yvar <- unlist(plot.structure[[plot_struct_idx]]$yvar)

    plot.list[[plot.list.idx]] <- single.plot(xvar, yvar, data.tibble)

    # add ylim
    if ("ylim" %in% names(plot.structure[[plot_struct_idx]]))
    {
        plot.list[[plot.list.idx]] <- plot.list[[plot.list.idx]] + ylim(
                unlist(
                        plot.structure[[plot.list.idx]]$ylim)
                )
    }
    
    plot.list[[plot.list.idx]] <- plot.list[[plot.list.idx]] + theme_classic()
    
    plot.list.idx <- plot.list.idx + 1
}

for (plot_struct_idx in 1:plot.structure.l)
{
    # get the (potential list of) y variable(s)
    # as this is a list and hence highly structured
    # hence, try to flatten it
    xvar <- unlist(plot.structure[[plot_struct_idx]]$xvar)
    yvar <- unlist(plot.structure[[plot_struct_idx]]$yvar)

    plot.list[[plot.list.idx]] <- single.plot(xvar, yvar, data.tibble.middle)

    # add ylim
    if ("ylim" %in% names(plot.structure[[plot_struct_idx]]))
    {
        plot.list[[plot.list.idx]] <- plot.list[[plot.list.idx]] + ylim(
                unlist(
                        plot.structure[[plot.list.idx]]$ylim)
                )
    }
    
    plot.list[[plot.list.idx]] <- plot.list[[plot.list.idx]] + theme_classic()
    
    plot.list.idx <- plot.list.idx + 1
}

base <- basename(file.name)
dir <- dirname(file.name)

distribution_file_name <- paste0(dir,"/",distribution_file_prefix,base)

sigmoidal_sex_det <- function(a, b, xseq) {
    return(data.frame(x = xseq,
            y = 1.0 / (1.0 + exp(-b * (xseq - a)))
            ))
}

sigmoidal_season <- function(a, b, tmax, xseq, tempseq) {
    return(data.frame(x = xseq,
            y = tmax / (1.0 + exp(-b * (tempseq - a)))
            ))
}

sigmoidal_plot_sex_det <- function(data, n.sample = 100) {

    # sample, say 100 sigmoidals from the data for 
    # sex determination
    n.sample <- 100
    rows.sampled <- sample(seq(1,nrow(data)),size = n.sample, replace=F)

    ampl.num <- as.numeric(params["amplitude"])
    xseq <- seq(-ampl.num,ampl.num, length.out = 100)
    
    list.sigmoidal.dfs <- list()

    for (i in 1:n.sample) {
        row <- data[rows.sampled[i],]
        sig_i <- sigmoidal_sex_det(a = row$a, 
                b = row$b,
                xseq)

        list.sigmoidal.dfs[[i]] <- sig_i
    }

    all.data <- bind_rows(list.sigmoidal.dfs,.id="id")

    p <- ggplot(data =all.data,
                mapping=aes(x=x, y=y,group=id)) +
                geom_line(alpha=0.5) +
                theme_classic()

    return(p)
} #sigmoidal_plot_sex_det

sigmoidal_plot_season <- function(data, n.sample = 100)
{
    # sample, say 100 sigmoidals from the data for 
    # sex determination
    n.sample <- 100
    rows.sampled <- sample(seq(1,nrow(data)),size = n.sample, replace=F)

    ampl.num <- as.numeric(params["amplitude"])
    
    max.t.season.num <- as.numeric(params["max_t_season"])

    # all the time steps over which the season fluctuates
    xseq <- seq(0, max.t.season.num, 
            length.out = 100)

    # now make the sinusoidal time series over this
    tempseq <- sin(xseq * 2 * pi / max.t.season.num)
    
    list.sigmoidal.dfs <- list()

    for (i in 1:n.sample) {
        row <- data[rows.sampled[i],]
        sig_i <- sigmoidal_season(a = row$temp_a, 
                b = row$temp_b,
                tmax = row$tmax,
                xseq, 
                tempseq)

        list.sigmoidal.dfs[[i]] <- sig_i
    }

    all.data <- bind_rows(list.sigmoidal.dfs,.id="id")

    envt.data <- data.frame(
            x = xseq,
            y = tempseq
            )

    p <- ggplot(data =all.data,
                mapping=aes(x=x, y=y)) +
                geom_line(mapping = aes(group = id), alpha=0.5) +
                geom_line(data = envt.data,
                        mapping = aes(x = x, y = y),
                        colour ="blue") +
                theme_classic()

    return(p)
}

title <- ""

if (exists("params") && "sf" %in% names(params))
{
    title <- paste0(
            "survival: ",params["sf"],", toptf: ",params["toptf"])
}


# then see whether there is a corresponding distribution file
if (file.exists(distribution_file_name))
{

    the_data_dist <- read.table(file=distribution_file_name, 
            sep=";",
            header=T)
    
    data_early <- the_data_dist[the_data_dist$time == min(the_data_dist$time),]
    data_late <- the_data_dist[the_data_dist$time == max(the_data_dist$time),]

    p_sex_det_early <- sigmoidal_plot_sex_det(data_early)
    p_sex_det_late <- sigmoidal_plot_sex_det(data_late)

    # put the sex det plot at the bottom
    # of the first row, which means it needs %
    plot.list <- append(
            plot.list, 
            list(p_sex_det_early,p_sex_det_late),
            after = floor(length(plot.list)/2)
            )

    p_season_early <- sigmoidal_plot_season(data_early)
    p_season_late <- sigmoidal_plot_season(data_late)

    # then put the season plot at the end
    plot.list[[length(plot.list) + 1]] <-  p_season_early
    plot.list[[length(plot.list) + 1]] <-  p_season_late
}

wrap_plots(plot.list,ncol=2, byrow=F) + plot_annotation(
        title=title)

file.name <- paste0("graph_",basename(file.name),".pdf")

ggsave(file.name,height= 3 * plot.structure.l, width = 20)


