
################################################################################
# Creates the initial simulated dataset in person-time format, assigning treatment
# and creating placeholders for time-varying variables including adherence and death.
# Time-varying variables will be subsequently updated iteratively with the 'create_sim_data' function.
# params: a list with the following structure,
## list(
##   sample_size = the number of individuals to simulate,
##   num_intervals = the number of intervals for each individual,
##   AKI_baseline = baseline risk of AKI,
##   AKI_adherent = absolute increase in risk of AKI if adherent to antihypertensive,
##   AKI_thiazide = absolute decrease in risk of AKI if adherent to thiazide,
##   abnormal_BP_baseline = baseline risk of abnormal blood pressure,
##   abnormal_BP_adherent = absolute decrease in risk of abnormal blood pressure if adherent to antihypertensive,
##   abnormal_BP_thiazide = further absolute decrease in risk of abnormal blood pressure if adherent to thiazide,
##   adherence_baseline = baseline probability of adherence,
##   adherence_thiazidecost = increase in probability of adherence due to cost of thiazides relative to ACEI,
##   adherence_AKI = decrease in probability of adherence due to AKI,
##   adherence_abnormal_BP = increase in probability of adherence due to abnormal blood pressure,
##   death_baseline = baseline risk of death,
##   death_abnormal_BP = increase in risk of death due to abnormal blood pressure,
##   death_adherence = decrease in risk of death due to adherence,
##   death_AKI = increase in risk of death due to AKI,
##   adherence_model = specify the causal model for the simulated data (1, 2, or 3)
##   )
################################################################################

initialize_sim_data <- function(params){
  
  with(params, {
    
    trt_ACEI <- as.numeric(c(rep(1, sample_size/2), rep(0, sample_size/2)))
    trt_thiazide <- 1-trt_ACEI
    
    out <- data.frame(
      id = rep(1:sample_size, each = num_intervals),
      trt_thiazide = rep(trt_thiazide, each = num_intervals),
      interval = rep(1:num_intervals, times = sample_size),
      AKI = rep(0, times = sample_size*num_intervals),
      abnormal_BP = rep(0, times = sample_size*num_intervals),
      adherent = rep(0, times = sample_size*num_intervals),
      dead = rep(0, times = sample_size*num_intervals)
    )
    
    return(out)
    
  })
  
}

################################################################################
# Initializes the simulated data then iteratively updates the time-varying variables
# params: a list with the following structure,
## list(
##   sample_size = the number of individuals to simulate,
##   num_intervals = the number of intervals for each individual,
##   AKI_baseline = baseline risk of AKI,
##   AKI_adherent = absolute increase in risk of AKI if adherent to antihypertensive,
##   AKI_thiazide = absolute decrease in risk of AKI if adherent to thiazide,
##   abnormal_BP_baseline = baseline risk of abnormal blood pressure,
##   abnormal_BP_adherent = absolute decrease in risk of abnormal blood pressure if adherent to antihypertensive,
##   abnormal_BP_thiazide = further absolute decrease in risk of abnormal blood pressure if adherent to thiazide,
##   adherence_baseline = baseline probability of adherence,
##   adherence_thiazidecost = increase in probability of adherence due to cost of thiazides relative to ACEI,
##   adherence_AKI = decrease in probability of adherence due to AKI,
##   adherence_abnormal_BP = increase in probability of adherence due to abnormal blood pressure,
##   death_baseline = baseline risk of death,
##   death_abnormal_BP = increase in risk of death due to abnormal blood pressure,
##   death_adherence = decrease in risk of death due to adherence,
##   death_AKI = increase in risk of death due to AKI,
##   adherence_model = specify the causal model for the simulated data (1, 2, or 3)
##   )
################################################################################

create_sim_data <- function(params) {
  
  with(params, {
    
    sim_data <- initialize_sim_data(params)
    
    sim_data$AKI[seq(1, 1+(num_intervals*(sample_size-1)), by=num_intervals)] <- 
      rbinom(n=sample_size, size=1,
             p = AKI_baseline)
    
    sim_data$abnormal_BP[seq(1, 1+(num_intervals*(sample_size-1)), by=num_intervals)] <- 
      rbinom(n=sample_size, size=1,
             p = abnormal_BP_baseline)
    
    sim_data$adherent[seq(1, 1+(num_intervals*(sample_size-1)), by=num_intervals)] <- 
      rbinom(n=sample_size, size=1,
             p = adherence_baseline)
    
    sim_data$dead[seq(1, 1+(num_intervals*(sample_size-1)), by=num_intervals)] <- 
      rbinom(n=sample_size, size=1,
             p = death_baseline +
               death_abnormal_BP*sim_data$abnormal_BP[seq(1, 1+(num_intervals*(sample_size-1)), by=num_intervals)] -
               death_adherence*sim_data$adherent[seq(1, 1+(num_intervals*(sample_size-1)), by=num_intervals)] +
               death_AKI*sim_data$AKI[seq(1, 1+(num_intervals*(sample_size-1)), by=num_intervals)])
    
    for (int in 2:num_intervals) {
      
      sim_data$AKI[seq(int, int+(num_intervals*(sample_size-1)), by=num_intervals)] <- 
        rbinom(n=sample_size, size=1,
               p = AKI_baseline + 
                 AKI_adherent*sim_data$adherent[seq((int-1), (int-1)+(num_intervals*(sample_size-1)), by=num_intervals)] -
                 AKI_thiazide*sim_data$adherent[seq((int-1), (int-1)+(num_intervals*(sample_size-1)), by=num_intervals)]*sim_data$trt_thiazide[seq((int-1), (int-1)+(num_intervals*(sample_size-1)), by=num_intervals)])
      
      sim_data$abnormal_BP[seq(int, int+(num_intervals*(sample_size-1)), by=num_intervals)] <- 
        rbinom(n=sample_size, size=1,
               p = abnormal_BP_baseline - 
                 abnormal_BP_adherent*sim_data$adherent[seq((int-1), (int-1)+(num_intervals*(sample_size-1)), by=num_intervals)] -
                 abnormal_BP_thiazide*sim_data$adherent[seq((int-1), (int-1)+(num_intervals*(sample_size-1)), by=num_intervals)]*sim_data$trt_thiazide[seq((int-1), (int-1)+(num_intervals*(sample_size-1)), by=num_intervals)])
      
      if (adherence_model == 1) {
        
        sim_data$adherent[seq(int, int+(num_intervals*(sample_size-1)), by=num_intervals)] <- 
          rbinom(n=sample_size, size=1,
                 p = adherence_baseline + 
                   adherence_thiazidecost*sim_data$adherent[seq((int-1), (int-1)+(num_intervals*(sample_size-1)), by=num_intervals)]*sim_data$trt_thiazide[seq((int-1), (int-1)+(num_intervals*(sample_size-1)), by=num_intervals)])
        
      }
      
      if (adherence_model == 2) {
        
        sim_data$adherent[seq(int, int+(num_intervals*(sample_size-1)), by=num_intervals)] <- 
          rbinom(n=sample_size, size=1,
                 p = adherence_baseline + 
                   adherence_thiazidecost*sim_data$adherent[seq((int-1), (int-1)+(num_intervals*(sample_size-1)), by=num_intervals)]*sim_data$trt_thiazide[seq((int-1), (int-1)+(num_intervals*(sample_size-1)), by=num_intervals)] - 
                   adherence_AKI*sim_data$AKI[seq(int, int+(num_intervals*(sample_size-1)), by=num_intervals)])
        
      }
      
      if (adherence_model == 3) {
        
        sim_data$adherent[seq(int, int+(num_intervals*(sample_size-1)), by=num_intervals)] <- 
        rbinom(n=sample_size, size=1,
               p = adherence_baseline + 
                 adherence_thiazidecost*sim_data$adherent[seq((int-1), (int-1)+(num_intervals*(sample_size-1)), by=num_intervals)]*sim_data$trt_thiazide[seq((int-1), (int-1)+(num_intervals*(sample_size-1)), by=num_intervals)] - 
                 adherence_AKI*sim_data$AKI[seq(int, int+(num_intervals*(sample_size-1)), by=num_intervals)] +
                 adherence_abnormal_BP*sim_data$abnormal_BP[seq(int, int+(num_intervals*(sample_size-1)), by=num_intervals)])
        
      }
      
      sim_data$dead[seq(int, int+(num_intervals*(sample_size-1)), by=num_intervals)] <- 
      rbinom(n=sample_size, size=1,
             p = death_baseline + 
               death_abnormal_BP*sim_data$abnormal_BP[seq(int, int+(num_intervals*(sample_size-1)), by=num_intervals)] -
               death_adherence*sim_data$adherent[seq(int, int+(num_intervals*(sample_size-1)), by=num_intervals)] +
               death_AKI*sim_data$AKI[seq(int, int+(num_intervals*(sample_size-1)), by=num_intervals)])
    
    }
    
    sim_data <- sim_data %>% 
      group_by(id) %>% 
      filter(cumsum(dead)==0 | (dead==1 & cumsum(dead)==1)) %>% 
      mutate(lag_adherent=lag(adherent,default=0),
             lag_abnormal_BP=lag(abnormal_BP,default=0),
             lag_AKI=lag(AKI,default=0)) %>%
      ungroup()
    
    return(sim_data)
    
  })
  
}