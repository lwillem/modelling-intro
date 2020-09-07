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

#' @title Plot of the population by health state
#'
#' @description  This function shows the spatial configuration of the population
#' with color codes for the health state and tracks one individual.
#'
#' @param pop_data        the vector with population data
#' @param area_size       the total area size
#' @param i_day           the current day
#' @param plot_time_delay the time delay between two plots
#'
#' @keywords external
#' @export
geo_plot_health_states <- function(pop_data,area_size,i_day,num_days,plot_time_delay)
{

  # if the time-delay is '0' ==>> skip figures
  if (plot_time_delay == 0){
    return(NULL)
  }

  # setup global variables for one participant 'X' (once!)
  if(!exists('participant_id')){

    # select all (centered) individuals (<1 from the centre)
    id_centered <- which(abs(pop_data$x_coord-(area_size/2)) < 1 &
                            abs(pop_data$y_coord-(area_size/2)) < 1)

    # sample one id and create matrix to log the x- and y-coordinates
    participant_id   <<- sample(id_centered,1)
    log_part_coord   <<- matrix(NA,nrow=2,ncol=num_days)
  }

  # (re)set figure layout on day 1
  if(i_day == 1){
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
       xlab = 'x coordinate (km)',
       ylab = 'y coordinate (km)',
       xlim = c(0,area_size),
       ylim = c(0,area_size+2),
       pch  = 2,
       main = paste('day',i_day));

  # add legend with color coding
  legend('topleft',
         c('S','I','R','V'),
         col  = 1:nlevels(pop_data_health_factor),
         pch  = 2,
         ncol = 4,
         cex  = legend_cex)

  # log coordinates of participant 'X' (adapt global variable)
  log_part_coord[,i_day]   <<- c(pop_data$x_coord[participant_id],pop_data$y_coord[participant_id])

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

  # add a day counter at the bottom left corner
  #text(area_size/2,area_size+1,paste('day',i_day),bg="white")
  #legend('top',paste('day',i_day),bg=NA, box.col = NA)

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
