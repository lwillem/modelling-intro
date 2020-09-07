############################################################################ #
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
############################################################################ #
#
# FUNCTION TO VISUALISE THE POPULATION IN THE RANDOM WALK TUTORIAL
#
############################################################################ #

#' @title Calculate the social contact probability
#'
#' @description  This function calculates the social contact probability based on
#' the average number of contacts per time step and the number of possible
#' social contacts at this time step.
#'
#' @note The maximum probability is limited to 0.999
#'
#' @param average_num_contacts   the average number of contacts per time step
#' @param num_possible_contacts  the number of possible contacts at this time step
#'
#' @keywords external
#' @export
get_contact_probability <- function(average_num_contacts,num_possible_contacts)
{

  # calculate the probability as the 'average' / 'possible'
  contact_probability <- average_num_contacts / num_possible_contacts

  # limit the probability to '0.999'
  if(contact_probability >= 1) {
    contact_probability <- 0.999
  }

  # return the probability
  return(contact_probability)

}


#' @title EXAMPLE to incorporate spatial vaccine refusal
#'
#' @description  This functions assumes spatial vaccine refusal in the
#' outer regions of the simulated area.
#'
#' @param pop_size          matrix with population data
#' @param vaccine_coverage  the vaccine coverage
#'
#' @keywords external
#' @export
sample_vaccine_refusal <- function(pop_data,vaccine_coverage){

  # (re)define the center of the simulated area
  area_size   <- max(c(pop_data$x_coord,pop_data$y_coord))
  area_center <- area_size / 2

  # (re)define population size
  pop_size <- nrow(pop_data)

  # define compliance radius
  radius <- (area_size/2) * vaccine_coverage * 0.9

  # select individuals in the central region, based on central x- and y-coordinates
  sel_x <- pop_data$x_coord < (area_center+radius) & pop_data$x_coord > (area_center-radius)
  sel_y <- pop_data$y_coord < (area_center+radius) & pop_data$y_coord > (area_center-radius)

  # combine the selection on x- and y-coordinate
  id_vaccine_potentials <- which(sel_x | sel_y)
  length(id_vaccine_potentials)

  # if we have to little vaccine potentials, add random individuals
  if(length(id_vaccine_potentials) < (pop_size*vaccine_coverage)){
    id_non_potentials        <- seq(1,pop_size) %in% id_vaccine_potentials
    required_potentials      <- (pop_size*vaccine_coverage) - length(id_vaccine_potentials)
    id_additional_potentials <- sample(id_non_potentials,required_potentials)
    id_vaccine_potentials    <- c(id_vaccine_potentials,id_additional_potentials)
  }

  # sample from the potential vaccineted individualss
  id_vaccinated <- sample(id_vaccine_potentials,pop_size*vaccine_coverage)


  # return indices
  return(id_vaccinated)
}



#' @title Print model parameters from the ibm_sirv_geo tutorial on the console
#'
#' @description  This functions provides an overview of the model parameters
#'
#' @note This function is created for the ibm_sirv_geo tutorial
#'
#' @keywords external
#' @export
print_sirv_geo_param <- function(){

  # collect possible parameter names
  all_param <- c('pop_size','num_days' ,'num_infected_seeds','vaccine_coverage',
                 'rng_seed','area_size','max_velocity','num_contacts_day',
                 'max_contact_distance', 'num_days_infected','transmission_prob',
                 'num_contacts_community_day','contact_prob_household','contact_prob_school',
                 'num_schools'
                 )

  # initiate string
  param_str <- ''

  # loop over the given parameter names, if present, add name & value
  for(i_param in all_param){
    if(exists(i_param)){
      param_str <- paste(param_str,'||',i_param,':',get(i_param))
    }
  }

  # print string
  cli::cat_rule('MODEL PARAMETERS')
  cli::cat_line(param_str,' ||')
  cli::cat_rule()
}

