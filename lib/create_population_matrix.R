#############################################################################
# This file is part of the SIMID course material
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright (C) 2019 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
#############################################################################
#
# FUNCTION TO CREATE A SYNTHETIC POPULATION
#
#############################################################################

#' @title Create a synthetic population with households
#'
#' @description This function creates a population with households
#'
#' @param pop_size  the final population size
#' @param area_size the size of the simulation area (to sample coordinates)
#'
#' @keywords external
#' @export
#pop_size <- 1e4
create_population_matrix <- function(pop_size,area_size)
{
  # demographic parameters
  ages_adult <- 18:60
  ages_child <- 1:18
  adult_age_tolerance     <- 0:5    # age tolerance between adults
  child_age_tolerance     <- 1:4    # age tolerance between children
  household_age_gap_min   <- 18     # min age gap between adults and children
  household_age_gap_max   <- 35     # max age gap age between adults and children

  # create the population
  pop_data         <- NULL  # start from empty matrix
  current_pop_size <- 0     # start from size '0'
  hh_id            <- 1     # a counter variable to track the household id

   # continue as long as 'population size' < 'target population size'
  while(current_pop_size<pop_size){

    # sample the age of adult 1
    age_adult1 <- sample(ages_adult, 1)

    # sample the age of adult 2, given adult 1
    age_adult2 <- sample(age_adult1 + adult_age_tolerance, 1)

    # get the possible child ages
    ages_child_option <- min(age_adult1,age_adult2) - (household_age_gap_min:household_age_gap_max )
    ages_child_option[!ages_child_option %in% ages_child]  <- NA
    ages_child_option <- c(NA,ages_child_option[!is.na(ages_child_option)])

    # sample the age of child 1
    age_child1 <- sample(ages_child_option, 1)

    # sample the age of child 2, given child 1
    age_child2 <- sample(age_child1 + child_age_tolerance, 1)

    # aggregate all ages with the household id
    hh_data <- data.frame(age = c(age_adult1,age_adult2,age_child1,age_child2),
                          hh_id = hh_id)

    # remove individuals with age 'NA' or negative ages (unborn)
    hh_data <- hh_data[!is.na(hh_data$age),]
    hh_data <- hh_data[hh_data$age>=0,]


    # add a household member id
    hh_data$member_id <- 1:nrow(hh_data)

    # add x- and y-coordindates
    hh_data$x_coord <- sample(seq(0,area_size,0.1),1) + sample(seq(-0.1,0.1,length=nrow(hh_data)))
    hh_data$y_coord <- sample(seq(0,area_size,0.1),1) + sample(seq(-0.1,0.1,length=nrow(hh_data)))

    # add hh_data to pop_data
    pop_data <- rbind(pop_data,
                      hh_data)

    # update statistics and household counter
    current_pop_size <- nrow(pop_data)
    hh_id    <- hh_id + 1

  } # end while-loop

  # select all individuals within the given population size
  pop_data <- pop_data[1:pop_size,]

  # inspect population characteristics
  head(pop_data)

  # add health state: susceptible
  pop_data <- data.frame(pop_data,
                         health              = 'S',           # column to store the health state
                         infector            = NA,            # column to store the source of infection
                         time_of_infection   = NA,            # column to store the time of infection
                         generation_interval = NA,            # column to store the generation interval
                         secondary_cases     = 0,             # column to store the number of secondary cases
                         stringsAsFactors = F)

  # create a figure with 6 subplots
  par(mfrow=c(2,3))
  hist(pop_data$age,-1:70,main='total population',xlab='age')
  hist(pop_data$age[pop_data$member_id==1],-1:70,main='adult 1',xlab='age')
  hist(pop_data$age[pop_data$member_id==2],-1:70,main='adult 2',xlab='age')
  hist(pop_data$age[pop_data$member_id==3],-1:70,main='child 1',xlab='age')
  hist(pop_data$age[pop_data$member_id==4],-1:70,main='child 2',xlab='age')
  hist(table(pop_data$hh_id),main='household size',xlab='household size')
  par(mfrow=c(1,1)) # restore the figure option: 1 plot per figure

  #geo_plot_health_states(pop_data,area_size,1,num_days,0.1)

  return(pop_data)

} # end function


#' @title Load a precomputed synthetic population with households
#'
#' @param pop_size  the final population size
#' @param area_size the size of the simulation area (to sample coordinates)
#'
#' @keywords external
#' @export
get_default_population_matrix <- function(pop_size,area_size)
{

  if(pop_size == 2e3 && area_size==20){

    print("Use pre-computed population")

    # use previously create population
    return(pop_matrix_2e3_20)

  } else{
    stop('The function "get_default_population_matrix" only works with pop_size = 2e3 and area_size=20')
  }

  return(NULL)
}


