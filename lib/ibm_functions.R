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
# Copyright (C) 2020 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
############################################################################ #
#
# FUNCTION TO VISUALISE THE POPULATION IN THE RANDOM WALK TUTORIAL
#
############################################################################ #

get_default_parameters <- function(){
  
  attach(list(pop_size = 1000,  # population size
              num_days = 50,    # time horizon
              num_infected_seeds = 3, # initial infections
              vaccine_coverage = 0.1, # vaccine state
              apply_spatial_vaccine_refusal  = TRUE, # are vaccine states randomly distributed
              
              area_size = 20,     # geo-spatial settings
              max_velocity = 0,
              
              num_days_infected  = 7, # disease parameter
              transmission_prob  = 0.1, # transmission dynamics
              num_contacts_day = 10, 
              max_contact_distance = 2,
              
              plot_time_delay  = 0.1,
              
              rng_seed = 2020))
  
  
}
#'
#' INDIVIDUAL-BASED MODEL (IBM) WITH:
#'   --> INDIVIDUAL-BASED WITH SPATIALLY EXPLICIT RANDOM WALK
#'   --> HEALTH STATES S, I, R, V
#'   --> VACCINE EFFICACY: 100%
#'   
run_ibm_random_walk <- function(pop_size = 1000,  # population size
                                num_days = 50,    # time horizon
                                num_infected_seeds = 3, # initial infections
                                vaccine_coverage = 0.1, # vaccine state
                                apply_spatial_vaccine_refusal  = FALSE, # is vaccination geographicaly clustered?
                                
                                area_size = 20,     # geo-spatial settings
                                max_velocity = 1,
                                
                                num_days_infected  = 7, # disease parameter
                                transmission_prob  = 0.1, # transmission dynamics
                                num_contacts_day = 10, 
                                max_contact_distance = 2,
                                
                                plot_time_delay  = 0, # visualisation parameter (0 = no plots)
                                
                                rng_seed = as.numeric(format(Sys.time(),'%S')) # random number seed = current time (seconds)
                                ){
  
  ######################################################### #
  # INITIALIZE POPULATION & MODEL PARAMETERS  ----
  ######################################################### #
  # initialize random number generator
  set.seed(rng_seed)
  
  # population vector: one row per individual, one column per attribute, row index = id
  pop_data     <- data.frame(health  = rep('S',length=pop_size),  # all individuals start in state 'S' (= susceptible)
                             x_coord = sample(seq(0,area_size,0.01),pop_size,replace = T), # sample random x coordinate
                             y_coord = sample(seq(0,area_size,0.01),pop_size,replace = T), # sample random y coordinate
                             infector            = NA,            # column to store the source of infection
                             time_of_infection   = NA,            # column to store the time of infection
                             generation_interval = NA,            # column to store the generation interval
                             secondary_cases     = 0,             # column to store the number of secondary cases
                             stringsAsFactors    = FALSE)         # option to treat characters as 'strings' instead of 'factors'
  
  # set vaccine coverage
  
  if(apply_spatial_vaccine_refusal){
    # option A: spatial clustering with respect to vaccine refusal
    id_vaccinated                  <- sample_vaccine_refusal(pop_data,vaccine_coverage)
  } else{
    # option B: random
    id_vaccinated                  <- sample(pop_size,pop_size*vaccine_coverage)
  }
  pop_data$health[id_vaccinated] <- 'V'
  
  # introduce infected individuals in the population
  id_infected_seeds                             <- sample(which(pop_data$health=='S'),num_infected_seeds)
  pop_data$health[id_infected_seeds]            <- 'I'
  pop_data$time_of_infection[id_infected_seeds] <- 0
  
  # set recovery parameters
  recovery_rate        <- 1/num_days_infected
  recovery_probability <- 1-exp(-recovery_rate)      # convert rate to probability
  
  # create matrix to log health states: one row per individual, one column per time step
  log_pop_data  <- matrix(NA,nrow=pop_size,ncol=num_days)
  
  # illustrate social contact radius
  if(plot_time_delay>0)
  geo_plot_social_contact_radius(pop_data,area_size,max_contact_distance,num_contacts_day,num_days)
  
  ####################################### #
  # RUN THE MODEL    ----                      
  ####################################### #
  
  # LOOP OVER ALL DAYS
  for(day_i in 1:num_days)
  {
    # step 1a: move at random [-1,1] units along the x and y axis
    step_vector      <- seq(-max_velocity,max_velocity,0.01)
    pop_data$x_coord <- pop_data$x_coord + sample(step_vector,pop_size,replace=T)
    pop_data$y_coord <- pop_data$y_coord + sample(step_vector,pop_size,replace=T)
    
    # step 1b: if an individual crossed the model world boundary: relocate at boundary
    pop_data$x_coord[pop_data$x_coord > area_size] <- area_size
    pop_data$y_coord[pop_data$y_coord > area_size] <- area_size
    pop_data$x_coord[pop_data$x_coord < 0]         <- 0
    pop_data$y_coord[pop_data$y_coord < 0]         <- 0
    
    # step 2: identify infected individuals
    boolean_infected <- pop_data$health == 'I'   # = boolean TRUE/FALSE
    ind_infected     <- which(boolean_infected)  # = indices
    num_infected     <- length(ind_infected)     # = number
    
    # step 3: calculate the distance matrix using the 'dist' function and store as matrix
    distance_matrix <- as.matrix(dist(pop_data[,c('x_coord','y_coord')],upper=T))
    
    # step 4: loop over all infected individuals
    p <- ind_infected[1]
    for(p in ind_infected)
    {
      # identify possible social contacts of person 'p'
      num_possible_contacts  <- sum(distance_matrix[p,] <= max_contact_distance)
      
      # calculate contact probability
      # tip: ?get_contact_probability
      contact_probability    <- get_contact_probability(num_contacts_day,num_possible_contacts)
      
      # new infections are possible if individuals are susceptible and within the range of the transmission distance
      flag_new_infection     <- pop_data$health == 'S' &
        distance_matrix[p,] <= max_contact_distance &
        rbinom(pop_size, size = 1, prob = contact_probability * transmission_prob)
      
      # mark new infected individuals
      pop_data$health[flag_new_infection] <- 'I'
      
      # log transmission details
      pop_data$infector[flag_new_infection]             <- p
      pop_data$time_of_infection[flag_new_infection]    <- day_i
      pop_data$secondary_cases[p]                       <- pop_data$secondary_cases[p] + sum(flag_new_infection)
      pop_data$generation_interval[flag_new_infection]  <- day_i - pop_data$time_of_infection[p]
    }
    
    # step 5: identify newly recovered individuals
    new_recovered <- boolean_infected & rbinom(pop_size, size = 1, prob = recovery_probability)
    pop_data$health[new_recovered] <- 'R'
    
    # step 6: log population health states
    log_pop_data[,day_i] <- pop_data$health
    
    # plot spatial configuration of the population by health state
    geo_plot_health_states(pop_data,area_size,day_i,num_days,plot_time_delay)
    
  } # end for-loop for each day
  
  ## PLOT RESULTS ----
  # reformat the log matrix with one row per individual and one column per time step
  # 'colSums' = sum per column
  log_s <- colSums(log_pop_data == 'S')  / pop_size
  log_i <- colSums(log_pop_data == 'I')  / pop_size
  log_r <- colSums(log_pop_data == 'R')  / pop_size
  log_v <- colSums(log_pop_data == 'V')  / pop_size
  
  
  # change figure configuration => 3 subplots
  #par(mfrow=c(1,3))
  
  m <- rbind(c(1, 1,1), c(2, 3, 4))
  layout(m)

  # final population overview
  geo_plot_health_states(pop_data,area_size,day_i,num_days,0.1,show_path = FALSE)
  
  # plot health states over time
  plot(log_s,
       type='l',
       xlab='Time (days)',
       ylab='Population fraction',
       main='Spatial IBM',
       ylim=c(0,1),
       lwd=2)
  lines(log_i,  col=2,lwd=2)
  lines(log_r,  col=3,lwd=2)
  lines(log_v,  col=4,lwd=2)
  
  legend('top',legend=c('S','I','R','V'),col=1:4,lwd=2,ncol=2,cex=0.7)
  
  boxplot(secondary_cases ~ time_of_infection, data=pop_data,
          xlab='time of infection (day)',
          ylab='secondary cases',
          main='secondary cases',
          ylim=c(0,10),
          xlim=c(0,num_days),
          xaxt='n')
  axis(1,seq(0,num_days,5))
  
  boxplot(generation_interval ~ time_of_infection, data=pop_data,
          xlab='time of infection (day)',
          ylab='generation interval (days)',
          main='generation interval',
          ylim=c(0,10),
          xlim=c(0,num_days),
          xaxt='n')
  axis(1,seq(0,num_days,5))

  ## PRINT MODEL PARAMETERS AND RESULTS ----
  # collect possible parameter names
  all_param <- c('pop_size','num_days' ,'num_infected_seeds','vaccine_coverage','apply_spatial_vaccine_refusal',
                 'rng_seed','area_size','max_velocity','num_contacts_day',
                 'max_contact_distance', 'num_days_infected','transmission_prob',
                 'num_contacts_community_day','contact_prob_household','contact_prob_school',
                 'num_schools','plot_time_delay'
  )
  
  print('MODEL PARAMETERS')
  
  # loop over the given parameter names, if present, add name & value
  for(i_param in all_param){
    if(exists(i_param)){
      # param_str <- paste(param_str,'||',i_param,':',get(i_param))
      print(paste(i_param,':',get(i_param)))
    }
  }
  
  # print total incidence
  print('-------------')
  print('MODEL RESULTS')
  
  print(paste0('total incidence: ',round((log_i[length(log_i)] + log_r[length(log_r)])*100,digits=2),'%'))
  
  # print peak details
  print(paste0('Peak prevalence: ',round(max(log_i)*100,digits=2),'%'))
  print(paste0('Peak day:        ',which(log_i == max(log_i)))[1])
}

#' @title Calculate the social contact probability
#'
#' @description  This function calculates the social contact probability based on
#' the average number of contacts per time step and the number of possible
#' social contacts at this time step.
#'
#' @note The maximum probability is limited to 0.95 (arbitrary choice)
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

  # limit the probability to '0.95'
  if(contact_probability >= 1) {
    contact_probability <- 0.95
  }

  # return the probability
  return(contact_probability)

}

run_ibm_location <- function(pop_size = 1000,  # population size
                              num_days = 50,    # time horizon
                              num_infected_seeds = 3, # initial infections
                              vaccine_coverage = 0.1, # vaccine state

                              num_days_infected  = 7, # disease parameter
                              transmission_prob  = 0.1, # transmission dynamics
                              num_contacts_day = 10, 
                              
                              plot_time_delay  = 0, # visualisation parameter (0 = no plots)
                              
                              rng_seed = as.numeric(format(Sys.time(),'%S')) # random number seed = current time (seconds)
){
  
  print("in progress...")
  
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


#' @title Plot of the population by health state
#'
#' @description  This function shows the spatial configuration of the population
#' with color codes for the health state and tracks one individual.
#'
#' @param pop_data        the vector with population data
#' @param area_size       the total area size
#' @param day_i           the current day
#' @param plot_time_delay the time delay between two plots
#'
#' @keywords external
#' @export
geo_plot_health_states <- function(pop_data,area_size,day_i,num_days,plot_time_delay,show_path=TRUE)
{
  
  # if the time-delay is '0' ==>> skip figures
  if (plot_time_delay == 0){
    return(NULL)
  }
  
  # (re)set figure layout on day 1
  if(day_i == 1){
    par(mfrow=c(1,1))
  }
  
  # clear the console
  flush.console()
  
  # set legend text size
  legend_cex <- 0.7
  
  # translate health states into a numeric order
  pop_data_health_factor <- factor(pop_data$health,levels=c('S','I','R','V'))
  
  # plot location and health state (color)
  plot(x    = pop_data$x_coord,
       y    = pop_data$y_coord,
       col  = pop_data_health_factor,
       xlab = 'x coordinate',
       ylab = 'y coordinate',
       xlim = c(0,area_size),
       ylim = c(0,area_size+2),
       pch  = 2,
       main = paste('day',day_i));
  
  # add legend with color coding
  legend('topleft',
         c('S','I','R','V'),
         col  = 1:nlevels(pop_data_health_factor),
         pch  = 2,
         ncol = 4,
         cex  = legend_cex)
  
  # track one individual? else ==>> skip 
  if (show_path){
    
    # setup global variables for one participant 'X' (once!)
    if(day_i == 1 | !exists('participant_id')){
      
      # select all (centered) individuals (<1 from the centre)
      id_centered <- which(abs(pop_data$x_coord-(area_size/2)) < 1 &
                             abs(pop_data$y_coord-(area_size/2)) < 1)
      
      # sample one id and create matrix to log the x- and y-coordinates
      participant_id   <<- sample(id_centered,1)
      log_part_coord   <<- matrix(NA,nrow=2,ncol=num_days)
    }
    
    # log coordinates of participant 'X' (adapt global variable)
    log_part_coord[,day_i]   <- c(pop_data$x_coord[participant_id],pop_data$y_coord[participant_id])
    log_part_coord <<- log_part_coord
    
    # add movement of participant 'X'
    lines(log_part_coord[1,],log_part_coord[2,],col=6,lwd=2)
    points(pop_data$x_coord[participant_id],
           pop_data$y_coord[participant_id],
           col=6,
           pch=2,
           lwd=5);
    
    # add legend for participant 'X'
    legend('topright',
           c('1 individual','path'),
           col  = 6,
           pch  = c(17,-1),
           lty  = c(0,1),
           ncol = 2,
           lwd  = 2,
           cex  = legend_cex)
  }
  
  # pause the system to make the time steps visual
  Sys.sleep(plot_time_delay)
  
} # end function


#' @title Plot the population and focus on the social contact radius
#'
#' @description This function shows the spatial configuration of the population
#' with color codes for the health state and the social contact radious of one individual.
#'
#' @param pop_data             the vector with population data
#' @param area_size            the total area size
#' @param max_contact_distance the max distance between 2 individuals for a contact
#' @param average_num_contacts the average number of contacts per day
#'
#' @keywords external
#' @export
geo_plot_social_contact_radius <- function(pop_data,area_size,max_contact_distance,average_num_contacts,num_days)
{
  
  # plot population
  geo_plot_health_states(pop_data,area_size,1,num_days,0.1)
  
  # add grid lines
  # note: 'abline' covers the full figure area and cannot be stoped at the model world boundary
  for(i_tick in 0:area_size){
    
    # define color and line type (i.e., solid for boundaries)
    line_col <- ifelse(i_tick %in% c(0,area_size),grey(0),grey(0.5))
    line_lty <- ifelse(i_tick %in% c(0,area_size),1,2)
    
    # plot horizontal and vertical lines
    lines(rep(i_tick,area_size+1),0:area_size,lty=line_lty,col=line_col)
    lines(0:area_size,rep(i_tick,area_size+1),lty=line_lty,col=line_col)
  }
  
  # get participant 'x'
  distance_matrix   <- as.matrix(dist(pop_data[,c('x_coord','y_coord')],upper=F,method = "euclidean"))
  distance_matrix[participant_id,participant_id] <- NA
  possible_contacts <- distance_matrix[participant_id,] < max_contact_distance
  points(pop_data$x_coord[possible_contacts],
         pop_data$y_coord[possible_contacts],
         pch=2,
         col="orange")
  
  # count possible contacts
  num_possible_contacts <- sum(possible_contacts,na.rm=T)
  
  # calculate the contact probability per possible contact
  contact_probability   <- get_contact_probability(average_num_contacts,num_possible_contacts)
  
  # set legend text size
  legend_cex <- 0.7
  
  legend('bottomleft',
         c(paste('average num. contacts:      ',average_num_contacts),
           paste('max. contact distance (km):',max_contact_distance),
           paste('possible num. contacts: ',num_possible_contacts),
           paste('contact probability:',trunc(contact_probability*100)/100)),
         col  = c(0,"orange",0),
         pch  = c(-1,2,-1),
         lty  = 0,
         ncol = 2,
         lwd  = 2,
         cex  = legend_cex)
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
  all_param <- c('pop_size','num_days' ,'num_infected_seeds','vaccine_coverage','is_vaccination_clustered',
                 'rng_seed','area_size','max_velocity','num_contacts_day',
                 'max_contact_distance', 'num_days_infected','transmission_prob',
                 'num_contacts_community_day','contact_prob_household','contact_prob_school',
                 'num_schools','plot_time_delay'
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
