#!/usr/bin/env python3
import pandas as pd
import os, re, sys

first = True
init = False
num_semicol_header = 0

# function to get the values of certain 
# columns at certain time points, meant to calculate 
# stuff pre and post perturbation
def fitness_at_time_points(file, maxrows, time_points, column_names):

    the_data = pd.read_csv(filepath_or_buffer=file,
                sep=";",
                nrows = maxrows)

    overall_data = None

    for time_point_i in time_points:

        subset = the_data.loc[the_data['time'] == time_point_i, column_names]

        subset.reset_index(drop=True, inplace=True)

        # append time points to column names
        colnames = list(subset.columns)
        for col_idx in range(0, len(colnames)):
            colnames[col_idx] = colnames[col_idx] + "_" + str(time_point_i)

        subset.columns = colnames

        overall_data = pd.concat([overall_data, subset], axis=1)

    return_str = overall_data.to_csv(path_or_buf=None,sep=";",index=False).strip()

    return_list_lines = return_str.split("\n")

    return(return_list_lines)
        


# analyze parameters at the end of the file
def analyze_parameters(lines,first=False):

    pars = {}

    for line in lines:
        mobj = line.split(";")

        if len(mobj) > 1:
            pars[mobj[0]] = mobj[1]

    return(pars)

def analyze_data(lines):

    data = [ [ float(celli) ] for celli in lines[0].split(";")[0:-1] ]

    # loop through the lines to collect the data
    for line in lines[1:]:
        splitline = line.split(";")[0:-1]

        for i in range(0,len(splitline)):
            data[i].append(float(splitline[i]))

    # now take averages

    avgs = []
    for i in range(0,len(data)):
        avgs.append(sum(data[i])/len(data[i]))

    return(avgs)

# processes the first line headers
# when making line headers for initial values
def process_first_line(line):

    # get the column names and split them into a list
    line_cols = line.strip().split(";")

    new_cols = ""

    for colname in line_cols:
        if not colname:
            continue

        new_cols += colname + "_t_0;" 

    return(new_cols)

def analyze_file(filename):

    global first;
    global num_semicol_header;

    global init;

    # indicator variable whether we are at first line
    # of the file to be read
    firstline = True

    flhead = ""

    lc = 0

    # indicator variable whether we are at the part
    # involving parameters
    parameter_part = False

    parameter_lines = []

    # the line where the parameter output
    # starts
    parline = -1

    # store the last line of data
    last_data_line = ""
    last_data_line_number = None

    # store the first line of data
    # in case we need initial values too
    first_data_line =""

    # the header of the resulting data file
    flhead = ""

    # open the file and read the stuff
    with open(filename) as infile:
        for line in infile:

            # see whether we also have to store the initial values
            if init:
                if firstline:
                    flhead += process_first_line(line)

            # update line count
            lc += 1

            # get the first line of data
            if lc == 2:
                first_data_line = line.strip()

            # if this is the first line store
            # the header
            if firstline:
                flhead += line.strip()
                firstline = False

            # if this is any other line starting
            # with a numerical value store the line
            # as it might be potentially the last one
            elif re.match(r"^\d",line):
                last_data_line = line
                last_data_line_number = lc

            # hold this until we have the parameter file
            if not parameter_part:
                if re.match("^\n",line) is not None:
                    parline = lc
                    parameter_part = True
                    parameter_lines += [line.strip()]
            elif parameter_part:
                parameter_lines += [line.strip()]

    if parline < 1:
        return

    parameters = analyze_parameters(parameter_lines)

    time_perturb = int(parameters["simulation_time_change"])
    season_time = int(parameters["max_t_season"])
    
    time_points = [time_perturb - 1, time_perturb + season_time - 1]

    # perform some extra operations in this case
    output_fitness = fitness_at_time_points(file = filename,
                                   maxrows = last_data_line_number - 1,
                                   time_points = time_points, 
                                   column_names = ["surviving_male_juvs", "surviving_female_juvs"])

    if len(output_fitness) < 2:
        return
        

    # prepare the data to be printed
    data = ""

    if init:

        # do error checking in terms of the csv values
#        count_semicol = len(re.findall(";", first_data_line))
#        count_semicol_data = len(re.findall(";", last_data_line))
#
#        assert(count_semicol == count_semicol_data)

        data += first_data_line.strip() 

    data += last_data_line.strip() 

    if first:

        header_line = ";".join(parameters.keys()) + ";" + flhead.strip() + output_fitness[0] + ";" + "file"

        # count number of occurrences of semicolon for error checking
        num_semicol_header = len(re.findall(";",header_line))

        print(header_line)

        first = False

    data_line = ";".join(parameters.values()) + ";" + data + output_fitness[1] + ";" + filename

    # count number of occurrences of semicolon for error checking
    num_semicol_data = len(re.findall(";",data_line))

    if (num_semicol_header != num_semicol_data):
        print(data_line)
        print(num_semicol_header)
        print(num_semicol_data)
        print("number of columns in data and header not equal...")
        sys.exit(1)

    print(data_line)


if len(sys.argv) > 2:

    # get initial values from the files as well
    init = True


# run the function on the indicated dir
for root, dir, files in os.walk(sys.argv[1]):

    for file in files:
        if re.match(r"(output|competition_output|sim|iter).*\d(.csv|.txt)?$",file) != None:
            analyze_file(os.path.join(root, file))
