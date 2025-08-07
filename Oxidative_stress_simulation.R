#'# ##' #'Questioning the utility of oxidative stress measurements as biomarkers of physiological condition and fitness
### Rachel R Reid, Davide M Dominoni, Jelle Boonekamp
##' 
##' 
##' Script aim:
#'Simulation aimed to test the consequences of low oxidative stress repeatability. 
#'11th November 2024
#'Updated: 7th August 2025


######Oxidative stress simulation 

# Load required libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(pwr)
library(mlmhelpr)
library(lme4)
library(plotly)
library(foreach)
library(doParallel)
library(reshape)



########
set.seed(123)

#####
# Simulation parameters for telomere dynamics and mortality


# Number of simulations to run per scenario
N=1000 


#Repeatability of oxidative stress across individuals
# Can be expanded to multiple levels for scenario testing 
repeatability_levels <- 0.15 

#Creates different scenarios of oxidative stress individual repeatability
#seq(0.15, 0.95, by = 0.35) 

# Causality strength levels for oxidative stress effect on telomere loss
causality_levels = c(0.6)

#Initial telomere parameters
avg_initial_telomere_length <- 30000 #Mean initial telomere length (bp)
sd_initial_telomere_length <- 8000 #Standard deviation around the initial telomere length adding individual variability
min_telomere_length = 1000 # min viable telomere length before mortality is triggered

#Telomere attrition due to cell division 
telomere_loss_by_cell_division <- 100 #Average annual telomere length (bp)
sd_telomere_loss_by_cell_division <- 20 #Variation around average annual telomere loss. 

# Oxidative stress scaling parameters
basestress<-1500 # Baseline oxidative stress value
sdstress<-450 # Standard deviation around oxidative stress
scale_varOS<-2 # Scaling factor for oxidative stress variability

# Mortality model parameters (Gompertz model)
R=0.02 # Baseline mortality rate (age-independent)
a=0.15# Age specific mortality rate (Gompertz slope)



#######
# Generating sample size scenarios
# Creates a vector of increasing sample sizes (quadratic growth)

i=1 # Initialize loop counter
ni=10 # Starting sample size
sample_size_levels<-ni # Initialize with first sample size
while ( ni < 860) {
  i<-i+1 # Increment loop counter
  ni<-ni+i*i # Increase sample size by i squared (quadratic growth)
  nii<-c(ni)
  sample_size_levels<-rbind(sample_size_levels,nii) # Amend new sample size to vector
}

# Number of repeated measures per individual
# This vector allows simulation across varying levels of sampling resolution 
nrep_levels<-c(0,1,2,3,4,9)

#####
# Helper Functions for simulating telomere dynamics and mortality across individuals

#####

# simulate_individual_with_mortality()
# Simulates life trajector of a single individual accounting for
# telomere dynamics, oxidative stress, and mortality risk using a Gompertz model
# Inputs:
# - repeatability: Individual variation in oxidative stress levels across years
# - causality: strength of oxidative stress effect on telomere loss
# Returns:
# - A dataframe containing yearly telomer length, oxidative stress, and survival status

#####

simulate_individual_with_mortality <- function(repeatability, causality) {
  
  # Initialize individual's parameters
  telomere_length <- rnorm(1, avg_initial_telomere_length, sd_initial_telomere_length) # Initial telomere length from normal distribution
  cumulative_ox_stress <- 0 # Tracks oxidative stress over time
  base_stress <-rnorm(1,basestress,sdstress) # Individual baseline oxidative stress
  ox_stress <- base_stress 
  #env$OS #Sets the oxidative stress level based on the individuals environment
  age <- 0 # Age starts at 0
  surv <- 1 #Survival status - (1 = alive)
  
  #Store yearly values in a dataframe
  df <- data.frame(age, telomere_length, cumulative_ox_stress, ox_stress) 
  ID<-as.factor(runif(1,0,1)) # Assign random unique ID
  
  
  # Simulate lifespan year-by-year until death
  while ( surv > 0) {
    age <- age + 1
    
    #variability <- sd_OS * scale_varOS*(1 - repeatability)  # Reduced variability with higher repeatability
    
    #ox_stress <- rnorm(1,(scale_meanOS+base_stress),variability) #Oxidative stress levels are drawn from a normal distribution and scaled using parameters assigned at the start of the script
    #Oxidative stress varies year to year depending on repeatability
    rn<-1-runif(1,repeatability,1)
    ox_stress<- base_stress+runif(1,-scale_varOS*basestress,scale_varOS*basestress)*(1-repeatability)
    
    # Telomere loss: stochastic cell division loss and oxidative stress related loss (scaled by causality)
    telomere_loss <- rnorm(1, telomere_loss_by_cell_division, sd_telomere_loss_by_cell_division) + ox_stress * causality 
    telomere_length <- max(0, telomere_length - telomere_loss)#Ensures telomere length never goes below 0
    cumulative_ox_stress <- cumulative_ox_stress + ox_stress
    
    # Apply Gompertz mortality or telomere based mortality
    if (runif(1) < R*exp(a*age) | telomere_length<min_telomere_length) {
      surv <- 0
    }
    
    df <- rbind(df, list(age, telomere_length, cumulative_ox_stress, ox_stress))
  }
  df<-cbind(ID,df)
  return(df)
}

#####
##Checking that the function works
#simulate_individual_with_mortality(0.3,1)

#####
# simulate_population_and_median_survival
# Simulates full population of individuals and computes:
# - Median survival time
# - Correlation between telomere length and oxidative stress at that time
# - Intra-class correlation coefficient (ICC) for oxidative stress
# - Statistical power to detect the telomere length and oxidative stress correlation
#
# Inputs
# - n: number of individuals in the population
# - nrep: number of repeated measures post-median survival
# - repeatability: within-individual variability in oxidative stress
# - causality: strength of oxidative stress effect on telomere loss
#Returns: 
# - Summary dataframe with correlation, ICC, power, and sample size
#####
simulate_population_and_median_survival <- function(n, nrep, repeatability, causality) {
  results <- replicate(
    n,
    simulate_individual_with_mortality(repeatability, causality),
    simplify = FALSE
  )
  
  # Combine individual data into a full population dataset
  d <- plyr::ldply(results, rbind) 
  
  # Calculate ICC for oxidative stress using a random intercept model
  fit <- lme4::lmer(ox_stress ~ (1 | ID), data = d)
  ICC <- icc(fit)[1, 2] 
  
  # Determine lifespan for each individual and calculate median survival time
  survival_times <- sapply(results, function(x) max(x$age)) 
  median_survival_time <- median(survival_times, na.rm = TRUE) 
  
  if (is.na(median_survival_time)) { 
    warning("Median survival time is NA.")
    return(list(NA, 0))
  }
  
  # Extract telomere and oxidative stress values for timepoints post median survival
  telomere_values <- melt(sapply(results, function(x) x$telomere_length[(median_survival_time):(median_survival_time+nrep)]))$value
  ID <- melt(sapply(results, function(x) x$ID[(median_survival_time):(median_survival_time+nrep)]))$value
  cumulative_ox_stress_values <- melt(sapply(results, function(x) x$cumulative_ox_stress[(median_survival_time):(median_survival_time+nrep)]))$value
  actual_ox_stress_values <- melt(sapply(results, function(x) x$ox_stress[(median_survival_time):(median_survival_time+nrep)]))$value
  dd<-data.frame(ID,telomere_values,actual_ox_stress_values)
  
  
#  correlation1 <- if (nrep == 0) {
#    cor(telomere_values, actual_ox_stress_values, use = "complete.obs")
#  } else {
#    #summary(lmer(scale(telomere_values)~ scale(actual_ox_stress_values)+(1|ID),data=dd))$ngrps
#    -sqrt(MuMIn::r.squaredGLMM(lmer(scale(telomere_values)~ scale(actual_ox_stress_values)+(1|ID),data=dd))[1])
#      }
  
  # Compute correlation between telomere length and oxidative stress
  correlation1 <-cor(telomere_values, actual_ox_stress_values, use = "complete.obs")
  
  # Estimate real sample size (accounting for repeats)
  real_N <- as.integer(dd%>%group_by(ID)%>%summarise(ID = first(ID)) %>% tally)
  
  # Estimate statistical power for the correlation at alpha = 0.05
  power <- pwr.r.test(r = correlation1, n = real_N, sig.level = 0.05)$power
  
  
  # Save individual-level data to CSV
  write.csv(d, paste0("PopulationData_Rep", repeatability, "_Caus", causality, ".csv"), row.names = FALSE)
  
  #Return population level summary stats
  df <- data.frame(correlation1, real_N, nrep, ICC, power) 
  return(df)
}

#####
# simulate_many_populations()
# Repeats the population level simulation multiple times (N)
# aggregating results to estimate average correlation, power, and ICC
#
#####
# Inputs:
# - N: number of repeated population simulations
# - n: individuals per population
# - nrep: number of repeated measurments per individual
# - repeatability: repeatability of oxidative stress
# - causality: strength of oxidative stress effect on telomere loss
# Returns:
# - Dataframe with average correlation, ICC, power, and standard errors. 
#####

simulate_many_populations <- function(N, n, nrep, repeatability,causality) {
  
  results <- replicate(
    N,
    simulate_population_and_median_survival(n, nrep, repeatability,causality),
    simplify = FALSE 
  )
  
  # Extract and summarise simulation outcomes across N replicates
  av_cor<- mean(sapply(results, function(x) mean(x$correlation1)))
  av_real_N<- mean(sapply(results, function(x) mean(x$real_N)))
  nrep<- mean(sapply(results, function(x) mean(x$nrep)))
  av_ICC<- mean(sapply(results, function(x) mean(x$ICC)))
  se_ICC <- sd(sapply(results, function(x) mean(x$ICC))) / sqrt(N)  # Compute SE for ICC
  av_power<- mean(sapply(results, function(x) mean(x$power)))
  se_power<-sd(sapply(results, function(x) mean(x$power)))/sqrt(N) # Computes SE for statistical power
  df<-data.frame(av_cor, av_ICC,se_ICC, av_real_N,nrep,av_power,se_power) #Returns dataframe
  return(df)
}

########
# Example simulation run
# Simulates 10 populations of 100 individuals each, with 0 repeat measurments, 
# no repeatability in oxidative stress, and full causality between oxidative stress and telomere loss
simulate_many_populations(N=10,n=100,nrep=0,repeatability=0,causality=1)

########
########

##########################
# Full simulation execution and visualisation 
##########################


#####
#Setup parallel processing
#Utilizes available CPU cores for faster execution of multiple parameter simulations
######

cl <- makeCluster(detectCores()-1) # Reserve one core for system processes
registerDoParallel(cl) # Register parallel backend

# Confirm parallel setup
getDoParWorkers()
getDoParName()



#####
# Create parameter grid
# Each row represents a unique combination of:
# - Repeatability of oxidative stress
# - Sample size per population
# - Causality of oxidative stress effect on telomere loss
# - Number of within-individual repeated samples
# - Number of population simulations per scenario
#####


# Define the parameter combinations
param_combinations <- expand.grid(
  Repeatability = repeatability_levels,
  Sample_Size = sample_size_levels,
  Causality = causality_levels,
  Nrep = nrep_levels,
  Ntot = N
)


#####
# Run simulations across all parameter combinations in parallel
####


start<- Sys.time()

power_simulation <- foreach(i = 1:nrow(param_combinations), .combine = rbind, .packages=c('lme4','mlmhelpr','dplyr','tidyr','pwr','reshape')) %dopar% {
  #Extract parameter set
  params <- param_combinations[i, ]
  
  #Run simulation for this parameter set
  data <- simulate_many_populations(N=params$Ntot, nrep = params$Nrep, n = params$Sample_Size, repeatability = params$Repeatability, causality = params$Causality)
  
  #Add metadata for plotting and tracking
  data$Repeatability <- params$Repeatability
  data$Sample_Size <- params$Sample_Size
  data$Nrep <- params$Nrep
  data$Causality <- params$Causality
  return(data)
}

print( Sys.time() - start)# Print total execution time


#####
# Post simulation data wrangling
#####

# Compute confidence intervals for power estimates (95% CI)
power_simulation<-power_simulation%>%mutate(ll=av_power-(qt(0.05,N)*se_power),ul=av_power+(qt(0.05,N)*se_power))

# Categorize ICC values into bins for visualization
power_simulation <- power_simulation %>%
  group_by(Nrep) %>%
  mutate(ICC_Bin = cut(av_ICC, breaks = c(-Inf, 0.3, 0.6, Inf))) %>%
  group_by(ICC_Bin) %>%
  mutate(Mean_ICC_Bin = mean(av_ICC, na.rm = TRUE)) %>% # Calculate mean ICC per bin
  ungroup() %>%
  mutate(ICC_Label = paste0("Mean ICC: ", round(Mean_ICC_Bin, 2))) # Format legend labels


# Compute 95% confidence intervals for ICC estimates
power_simulation <- power_simulation %>%
  mutate(
    ICC_lower = av_ICC - qt(0.975, N-1) * se_ICC,
    ICC_upper = av_ICC + qt(0.975, N-1) * se_ICC
  )



#####
# Plotting results
#####

# Set consistent publication-quality plot theme
jelle<- theme_update(
  panel.grid.major=element_line(colour=NA),
  panel.grid.minor=element_line(colour=NA),
  panel.background = element_rect(colour =NA ,fill=NA,linewidth=0.5),
  axis.title.x=element_text(size=20,hjust=0.5,vjust=0.5,angle=0),
  axis.title.y=element_text(size=20,hjust=0.5,vjust=0.5,angle=90),
  axis.text.x=element_text(colour="black",angle=0,size=10),
  axis.text.y=element_text(colour="black",angle=0,size=10),
  axis.ticks=element_line(colour="black",linewidt=0.5),
  axis.line.x=element_line(linewidth=0.5),
  axis.line.y=element_line(linewidth=0.5))


######
# Plot 1: Power vs Sample size, grouped by ICC levels
#####

ggplot(power_simulation, aes(x = av_real_N, y = av_power, color = ICC_Label, fill = ICC_Label)) +
  geom_line(size = 1) +  # Power curve
  geom_ribbon(aes(ymin = ll, ymax = ul), alpha = 0.2, linetype = "blank") +  # Shaded confidence intervals
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +  # Power threshold
  labs(
    title = "Sample Size vs Power for Detecting Correlation at Median Survival Time",
    x = "Sample Size",
    y = "Power",
    color = "Mean ICC per Bin",
    fill = "Mean ICC per Bin") +
  scale_x_continuous(limits=c(0, 450), breaks=c(seq(0, 450, by=50)))+
  scale_y_continuous(limits=c(0, 1), breaks=c(seq(0, 1, by=0.25)))+
  facet_wrap(~Nrep+Causality, ncol=3,labeller = 
               labeller(
                 Nrep = ~ paste("n repeated samples: ", .),
                 Causality = ~ paste("causality: ", .),
                 .multi_line = FALSE
               ))


#####
# Plot 2: Power curves under low repeatability (ICC <0.2), by nrep
#####

ss<-subset(power_simulation, power_simulation$Mean_ICC_Bin<0.2)

ggplot(ss, aes(x = av_real_N, y = av_power, color = as.factor(Nrep+1), fill = as.factor(Nrep+1))) +
  geom_line(size = 1) +  # Power curve
  geom_ribbon(aes(ymin = ll, ymax = ul), alpha = 0.2, linetype = "blank") +  # Shaded confidence intervals
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +  # Power threshold
  labs(
    title = "Low repeatability scenario (ICC = 0.2); 60% causality",
    x = "Sample Size",
    y = "Power",
    color = "number of within-individual samples",
    fill = "number of within-individual samples") +
  scale_x_continuous(limits=c(50, 350), breaks=c(seq(0, 450, by=50)))+
  scale_y_continuous(limits=c(0.25, 1), breaks=c(seq(0, 1, by=0.25)))+
  theme_minimal()


#####
# Plot 3: Relationship between repeated measures and sample sizse needed to reach 80% power
####

sss<-subset(power_simulation, power_simulation$Causality==0.6 & power_simulation$av_power<0.85 &power_simulation$av_power>0.79)%>%group_by(nrep)%>%summarise(mean_av_real_N=mean(av_real_N))


ggplot(sss, aes(x = nrep+1, y = mean_av_real_N)) +
  geom_point(size = 5) +  # Power curve
  geom_smooth(method=lm, formula=(y~x*exp(-(x))))+
  labs(
    title = "Sample Size vs repeated measures for achieving power >80%",
    x = "Number of repeated samples",
    y = "Sample Size") +
  scale_x_continuous(limits=c(1, 11), breaks=c(seq(0, 11, by=1)))+
  scale_y_continuous(limits=c(170, 300), breaks=c(seq(170, 300, by=30)))












