require('rmarkdown')

project_dir <- path.expand(".")
Rmd_dir <- file.path(project_dir, "Rmd")
html_dir <- file.path(project_dir, "docs")

Rmd_files <- c("index","intro_to_r", "intro_to_SIR",
               "intro_to_RStudio",'intro_to_IBM_coding',
               'intro_to_IBM_walk','intro_to_IBM_location',
               'HTA_flu_basic','HTA_flu_CEAresults',
               'HTA_flu_uncertainty','HTA_flu_evppi') 

for (Rmd_file in Rmd_files) {
  render(file.path(Rmd_dir, sprintf("%s.Rmd", Rmd_file)),
         output_dir = html_dir)
}



# set WTP range, and step size
opt_wtp <- seq(0,40000,2000)

# initialise matrix 
inmb_wtp_summary <- data.frame(wtp = opt_wtp,
                               inmb = NA,
                               ceac = NA,
                               ceaf = NA,
                               inmb_perfect_info = NA)

i <- 2
inmb_wtp_summary[i,]
for(i in 1:length(opt_wtp)){
  # calculate the inmb per sample
  inmb_sample <- (QALY_gained_disc * opt_wtp[i]) - incr_total_cost
  
  # get the expected INMB 
  inmb_wtp_summary$inmb[i] <- mean(inmb_sample)
  
  # inspect the number of samples with a positive INMB
  inmb_wtp_summary$ceac[i] <- sum(inmb_sample>=0) / nsamples
  
  # INMB with perfect information (i.e. set negative to 0)
  if(inmb_wtp_summary$inmb[i] >= 0){
    inmb_sample[inmb_sample < 0] <- 0
    inmb_wtp_summary$inmb_perfect_info[i] <- mean(inmb_sample)
  } else{
    inmb_sample[inmb_sample > 0] <- 0
    inmb_wtp_summary$inmb_perfect_info[i] <- mean(inmb_sample)
  }
}

# CEAC
plot(x = inmb_wtp_summary$wtp,
     y = inmb_wtp_summary$ceac,
     xlab = 'Willingness To Pay',
     ylab = 'Probability cost effective',
     main = 'Cost-Effectiveness Acceptibility Curve',
     ylim = 0:1)

# CEAF
inmb_wtp_summary$ceaf[inmb_wtp_summary$inmb > 0] = inmb_wtp_summary$ceac[inmb_wtp_summary$inmb > 0]
plot(x = inmb_wtp_summary$wtp,
     y = inmb_wtp_summary$ceaf,
     xlab = 'Willingness to pay',
     ylab = 'Probability to be cost-effective',
     main = 'Cost-Effectiveness Acceptibility Frontier',
     ylim = 0:1)

# Combination CEAC & CEAF
plot(x = inmb_wtp_summary$wtp,
     y = inmb_wtp_summary$ceac,
     type = 'l',
     lty = 3,
     xlab = 'Willingness To Pay',
     ylab = 'Probability to be cost-effective',
     main = 'Cost-Effectiveness Acceptibility',
     ylim = 0:1)
lines(x = inmb_wtp_summary$wtp,
      y = inmb_wtp_summary$ceaf,
      type = 'l',
      lty = 1,
      lwd = 3)
legend('topleft',
       c('CEAC', 'CEAF'),
       lty = c(3,1),
       lwd = c(1,3))

# NET-LOSS
inmb_wtp_summary$net_loss <- abs(inmb_wtp_summary$inmb_perfect_info - inmb_wtp_summary$inmb)

plot(x = inmb_wtp_summary$wtp,
     y = inmb_wtp_summary$net_loss,
     type = 'l',
     xlab = 'Willingness To Pay',
     ylab = 'Expected net-loss',
     main = "Expected Value of Perfect Information")

#################### #

set.seed(131)
n_samples <- 1000

# cost
c <- vector(mode = "list", length = 6)
names(c) <- c("Strategy 1, Grp 1", "Strategy 1, Grp 2", "Strategy 2, Grp 1",
              "Strategy 2, Grp 2", "Strategy 3, Grp 1", "Strategy 3, Grp 2")
c[[1]] <- rlnorm(n_samples, 2, .1)
c[[2]] <- rlnorm(n_samples, 2, .1)
c[[3]] <- rlnorm(n_samples, 11, .15)
c[[4]] <- rlnorm(n_samples, 11, .15)
c[[5]] <- rlnorm(n_samples, 11, .15)
c[[6]] <- rlnorm(n_samples, 11, .15)

# effectiveness
e <- c
e[[1]] <- rnorm(n_samples, 8, .2)
e[[2]] <- rnorm(n_samples, 8, .2)
e[[3]] <- rnorm(n_samples, 10, .8)
e[[4]] <- rnorm(n_samples, 10.5, .8)
e[[5]] <- rnorm(n_samples, 8.5, .6)
e[[6]] <- rnorm(n_samples, 11, .6)

# cost and effectiveness by strategy and simulation
library("data.table")
ce <- data.table(sample = rep(seq(n_samples), length(e)),
                 strategy = rep(paste0("Strategy ", seq(1, 3)), 
                                each = n_samples * 2),
                 grp = rep(rep(c("Group 1", "Group 2"),
                               each = n_samples), 3),
                 cost = do.call("c", c), qalys = do.call("c", e))
head(ce)

library("hesim")
ktop <- 200000
cea_out <-  cea(ce, k = seq(0, ktop, 500), sample = "sample", strategy = "strategy",
                grp = "grp", e = "qalys", c = "cost")

cea_pw_out <-  cea_pw(ce,  k = seq(0, ktop, 500), comparator = "Strategy 1",
                      sample = "sample", strategy = "strategy", grp = "grp",
                      e = "qalys", c = "cost")

library("magrittr") # Use pipes
icer(cea_pw_out, k = 50000) %>%
  format()

head(cea_pw_out$delta)

library("ggplot2")
theme_set(theme_bw())
plot_ceplane(cea_pw_out, k = 50000)
plot_ceac(cea_out)
plot_ceac(cea_pw_out)
plot_ceaf(cea_out)
plot_evpi(cea_out)


n_samples <- length(comparator_total_cost)
ktop <- 5e4
# cost and effectiveness by strategy and simulation
library("data.table")
ce <- data.table(sample = rep(seq(n_samples), length(e)),
                 strategy = rep(paste0("Strategy ", seq(1, 2)), 
                                each = n_samples),
                 grp = rep(rep(c("Group 1"),
                               each = n_samples), 2),
                 cost = c(comparator_total_cost,
                          program_total_cost), 
                 qalys = c(-comparator_QALYloss_disc,
                            -program_QALYloss_disc))
head(ce)


cea_pw_out <-  cea_pw(ce,  k = seq(0, ktop, 500), comparator = "Strategy 1",
                      sample = "sample", strategy = "strategy", grp = "grp",
                      e = "qalys", c = "cost")
library("magrittr") # Use pipes
icer(cea_pw_out, k = 50000) %>%
  format()

plot_ceplane(cea_pw_out, k = 25000)
plot_ceac(cea_pw_out)

ce_incr <- data.table(sample = rep(seq(n_samples), length(e)),
                 strategy = rep(paste0("Strategy ", seq(1, 2)), 
                                each = n_samples),
                 grp = rep(rep(c("Group 1"),
                               each = n_samples), 2),
                 cost = c(0*comparator_total_cost,
                          incr_total_cost), 
                 qalys = c(0*comparator_QALYloss_disc,
                           QALY_gained_disc))

cea_out <-  cea(ce_incr,  k = seq(0, ktop, 500),
                sample = "sample", strategy = "strategy", grp = "grp",
                e = "qalys", c = "cost")

plot_ceac(cea_out)
plot_ceaf(cea_out)
plot_evpi(cea_out)
?plot_ceaf

plot(comparator_total_cost,comparator_QALYloss_disc, xlim = range(0,program_total_cost))
points(program_total_cost,program_QALYloss_disc, col = 2)
icer(cea_pw_out)

mean(program_total_cost - comparator_total_cost) / 
mean(comparator_QALYloss_disc - program_QALYloss_disc)

