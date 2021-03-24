#!/usr/bin/env Rscript 

#--vanilla

library("ggplot2")
library("gridExtra")
library("tidyr")

# get command line arguments
args = commandArgs(trailingOnly=TRUE)

# give an error message if you do not provide it with a simulation file name
if (length(args) < 1)
{
    print("provide a simulation file name")
    stop()
}


# find out data line
find_out_data_line <- function(filename) {

    f <- readLines(filename)

    seq_lines <- seq(1,length(f),1)

    for (line_i in seq_lines)
    {
        if (length(grep("^generation",f[[line_i]])) > 0)
        {
            return(line_i)
        }
    } # end for
}

# find out where the parameter listing starts
# so that we can read in the data part of the file 
# without having it messed up by the subsequent parameter listing
find_out_param_line <- function(filename) {

    f <- readLines(filename)

    # make a reverse sequence
    seqq <- seq(length(f),1,-1)

    # go through each line in the data file and find first line
    # where data is printed (i.e., a line which starts with a digit)
    for (line_i in seqq)
    {
        print(f[[line_i]])
        print(line_i)
        if (length(grep("^\\d",f[[line_i]])) > 0)
        {
            return(line_i)
        }
    }

    return(NA)
}

# get the line at which the data starts
data_row <- find_out_data_line(args[1])

if (is.na(data_row))
{
    print("cannot find data...")
    stop()
}

# read in data frame of corresponding simulation
the.data <- read.table(args[1], header=T, skip = data_row - 1, sep=";")

# now use ggplot2 to plot stuff
# put the.data in a format so that trait combines
# b, sr1 and sr2
the.data.t <- pivot_longer(the.data
        ,cols=c("b","sr1","sr2")
        ,names_to="trait"
        ,values_to="trait_values")


p1.a <- ggplot(data=the.data.t,
            aes(x=generation, y=trait_values, colour=trait)) +
            geom_line() +
            theme_classic() + 
            labs(x =""
                    ,y="Trait values")

### Dispersal ####

the.data.t <- pivot_longer(the.data
        ,cols=c("df","dm")
        ,names_to="trait"
        ,values_to="trait_values")

p1.b <- ggplot(data=the.data.t,
            aes(x=generation, y=trait_values, colour=trait)) +
            geom_line() +
            theme_classic() + 
            labs(x =""
                    ,y="Dispersal")

### Consanguinity ####
names.q <- grep("^Q",names(the.data),value=T)

the.data.t <- pivot_longer(the.data
        ,cols=all_of(names.q)
        ,names_to="Consanguinity"
        ,values_to="trait_values")

p2 <- ggplot(data=the.data.t
        ,aes(x=generation, y=trait_values, colour=Consanguinity)) +
            geom_line() +
            theme_classic() + 
            labs(x =""
                    ,y="Q")

            
### Reproductive values ####
names.q <- grep("^v",names(the.data),value=T)
print(names.q)

the.data.t <- pivot_longer(the.data
        ,cols=all_of(names.q)
        ,names_to="RV"
        ,values_to="trait_values")

p3 <- ggplot(data=the.data.t
        ,aes(x=generation, y=trait_values, colour=RV)) +
            geom_line() +
            theme_classic() + 
            labs(x =""
                    ,y="Reproductive value")

### Class frequencies ####
names.q <- grep("^u",names(the.data),value=T)
print(names.q)

the.data.t <- pivot_longer(the.data,
        ,cols=all_of(names.q)
        ,names_to="u"
        ,values_to="trait_values")

p4 <- ggplot(data=the.data.t
        ,aes(x=generation, y=trait_values, colour=u)) +
            geom_line() +
            theme_classic() + 
            labs(x =""
                    ,y="Class frequencies")

### lambda ####
p5 <- ggplot(data=the.data
        ,aes(x=generation, y=lambda)) +
            geom_line() +
            theme_classic() + 
            labs(x =""
                    ,y="Dominant eigenvalue")

big_plot <- arrangeGrob(p1.a, p1.b, p2, p3, p4, p5, nrow=6,ncol=1)
the.base.name <- basename(args[1])

output_file_name <- paste(
        "graph_"
        ,the.base.name
        ,".pdf"
        ,sep="")

ggsave(output_file_name, big_plot, height = 20)

