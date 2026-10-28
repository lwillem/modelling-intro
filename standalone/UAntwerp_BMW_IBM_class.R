# This file is part of the course material of "Health economics 
# and mathematical models of infectious diseases (2018FBDBMW)" 
# as part of the Master of Biomedical Sciences. 
#
# For more info: see https://lwillem.github.io/modelling-intro
# 
# Author: Prof. Lander Willem (2026)
####################################################################### #

## SETUP ----

# install package to download the help function
install.packages('curl')

# load package 'curl'
library('curl')

# download R file with IBM function to your computer (in your working directory)
curl_download("https://github.com/lwillem/modelling-intro/raw/master/lib/ibm_functions.R", destfile = "ibm_functions.R")

# load the IBM functions
source('ibm_functions.R')

## RANDOM WALK ----

# run the core function with default parameters
# note: the model parameters that are printed to the console, can be used as function argument to change a setting
run_ibm_random_walk()

# run IBM with random walk (see console)
run_ibm_random_walk(plot_time_delay = 0.2,
                    rng_seed = 2000)

# run IBM with fixed population (see console)
run_ibm_random_walk(plot_time_delay = 0.2,
                    max_step_size = 0,
                    rng_seed = 2000)

## RANDOM WALK: scenarios ----
# 1. number of infected seeds to start the epidemic
run_ibm_random_walk(num_infected_seeds = 10, # change number of infected seeds to 10 (instead of 3)
                    rng_seed = 2000,         # fix the random number generator stream
                    add_baseline = TRUE)     # to include a run with default settings

# 2. population size
run_ibm_random_walk(pop_size = 2000,     # change population size to 100 (instead of 1000)
                    rng_seed = 2000,     # fix the random number generator stream
                    add_baseline = TRUE) # to include a run with default settings

# 3. no movement
run_ibm_random_walk(max_step_size = 0,   # change maximum step-size to 5 (instead of 2)
                    rng_seed = 2000,     # fix the random number generator stream
                    add_baseline = TRUE) # to include a run with default settings

# 3. a lot of movement
run_ibm_random_walk(max_step_size = 5,   # change maximum step-size to 5 (instead of 2)
                    rng_seed = 2000,     # fix the random number generator stream
                    add_baseline = TRUE) # to include a run with default settings
# caution: model artefact!

## LOCATION SPECIFIC MIXING ----
# see https://lwillem.github.io/modelling-intro/IBM_location.html

# inspect household set-up (by omitting the brackets, no function call)
head(create_population_matrix, 9)

# inspect the set_schools function (by omitting the brackets, no function call)
set_schools

# inspect the set_workplaces function (by omitting the brackets, no function call)
set_workplaces

# run the location-specific IBM with default parameters (see console)
run_ibm_location()

## LOCATION SPECIFIC MIXING: scenarios ----
# 1. number of infected seeds to start the epidemic
run_ibm_location(num_infected_seeds = 10,        # change number of infected seeds to 10 (instead of 3)
                 bool_show_demographics = FALSE, # do not show the demographics plot
                 rng_seed = 2000,                # fix the random number generator stream
                 add_baseline = TRUE)            # to include a run with default settings

# 2. population size
run_ibm_location(pop_size = 4000,                # change population size to 4000 (instead of 2000)
                 bool_show_demographics = FALSE, # do not show the demographics plot
                 rng_seed = 2000,                # fix the random number generator stream
                 add_baseline = TRUE)            # to include a run with default settings

# 3. no school contacts
run_ibm_location(contact_prob_school = 0,        # change contact behaviour
                 bool_show_demographics = FALSE, # do not show the demographics plot
                 rng_seed = 2000,                # fix the random number generator stream
                 add_baseline = TRUE)            # to include a run with default settings

# 4. no workplace contacts
run_ibm_location(contact_prob_work = 0,          # change contact behaviour
                 bool_show_demographics = FALSE, # do not show the demographics plot
                 rng_seed = 2000,                # fix the random number generator stream
                 add_baseline = TRUE)            # to include a run with default settings

# 4. no household contacts
run_ibm_location(contact_prob_household = 0,     # change contact behaviour
                 bool_show_demographics = FALSE, # do not show the demographics plot
                 rng_seed = 2000,                # fix the random number generator stream
                 add_baseline = TRUE)            # to include a run with default settings
