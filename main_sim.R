source('./src_sim/dependencies.R')
source('./src_sim/functions.R')


######################################################
############## Generate simulated data  ##############
######################################################

#Set sim parameters
params <- list(
  sample_size = 1000000,
  num_intervals = 24,
  AKI_baseline = 0.05,
  AKI_adherent = 0.035,
  AKI_thiazide = 0.035,
  abnormal_BP_baseline = 0.95,
  abnormal_BP_adherent = 0.6,
  abnormal_BP_thiazide = 0.15,
  adherence_baseline = 0.6,
  adherence_thiazidecost = 0.2,
  adherence_AKI = 0.5,
  adherence_abnormal_BP = 0.2,
  death_baseline = 0.035,
  death_abnormal_BP = 0.01,
  death_adherence = 0.03,
  death_AKI = 0.01
)

sim_data <- list(
  create_sim_data(append(params, list(adherence_model=1, separate_ZA_ZY=F))),
  create_sim_data(append(params, list(adherence_model=2, separate_ZA_ZY=F))),
  create_sim_data(append(params, list(adherence_model=3, separate_ZA_ZY=F))),
  create_sim_data(append(params, list(adherence_model=1, separate_ZA_ZY=T))),
  create_sim_data(append(params, list(adherence_model=2, separate_ZA_ZY=T))),
  create_sim_data(append(params, list(adherence_model=3, separate_ZA_ZY=T)))
)
  
save(sim_data, file='./sim_data.Rdata')

#######################################################
# Compute weights for intervention on prior death
# which can be used for covariate balance comparisons 
# However, we do not use these in the manuscript 
# Instead we present covariate and adherence comparisons
# conditional on past survival
#######################################################

for (d in 1:length(sim_data)) {
  
  sim_data[[d]]$Wt_death <- 1-predict(lm(dead ~ abnormal_BP + adherent + AKI, data=sim_data[[d]]), 
                                      newdata=sim_data[[d]])
  
  sim_data[[d]] <- sim_data[[d]] %>% group_by(id) %>%
    mutate(Wt_death = cumprod(1/lag(Wt_death, default=1))) %>%
    ungroup()
  
}

######################################################
################ Total (ITT) effect  #################
######################################################

total_dead_1 <- sim_data[[1]] %>% group_by(interval, trt_thiazide) %>% summarise(cond_risk=mean(dead)) %>% 
  group_by(trt_thiazide) %>% mutate(CI=1-cumprod(1-cond_risk)) %>% 
  {bind_rows(as.data.frame(.),data.frame(interval=c(0,0),trt_thiazide=c(0,1),cond_risk=c(0,0),CI=c(0,0)))} %>%
  ggplot(data=., aes(x=interval, y=CI, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(0,24)) + ylim(c(0,0.5)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('Cumulative incidence') +
  xlab('') + theme(text = element_text(size = 10))

####################################################################################################

total_dead_2 <- sim_data[[2]] %>% group_by(interval, trt_thiazide) %>% summarise(cond_risk=mean(dead)) %>% 
  group_by(trt_thiazide) %>% mutate(CI=1-cumprod(1-cond_risk)) %>% 
  {bind_rows(as.data.frame(.),data.frame(interval=c(0,0),trt_thiazide=c(0,1),cond_risk=c(0,0),CI=c(0,0)))} %>%
  ggplot(data=., aes(x=interval, y=CI, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(0,24)) + ylim(c(0,0.5)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('Cumulative incidence') +
  xlab('') + theme(text = element_text(size = 10))

####################################################################################################

total_dead_3 <- sim_data[[3]] %>% group_by(interval, trt_thiazide) %>% summarise(cond_risk=mean(dead)) %>% 
  group_by(trt_thiazide) %>% mutate(CI=1-cumprod(1-cond_risk)) %>% 
  {bind_rows(as.data.frame(.),data.frame(interval=c(0,0),trt_thiazide=c(0,1),cond_risk=c(0,0),CI=c(0,0)))} %>%
  ggplot(data=., aes(x=interval, y=CI, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(0,24)) + ylim(c(0,0.5)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('Cumulative incidence') +
  xlab('') + theme(text = element_text(size = 10))

######################################################
######### Adherence balance for total effect  ########
######################################################

total_adherence_1 <- sim_data[[1]] %>% group_by(interval, trt_thiazide) %>% summarise(adherence=mean(adherent)) %>% 
  ggplot(data=., aes(x=interval, y=adherence, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(1,24)) + ylim(c(0,1)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('Proportion') +
  xlab('') + theme(text = element_text(size = 10))

####################################################################################################

total_adherence_2 <- sim_data[[2]] %>% group_by(interval, trt_thiazide) %>% summarise(adherence=mean(adherent)) %>% 
  ggplot(data=., aes(x=interval, y=adherence, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(1,24)) + ylim(c(0,1)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('Proportion') +
  xlab('') + theme(text = element_text(size = 10))

####################################################################################################

total_adherence_3 <- sim_data[[3]] %>% group_by(interval, trt_thiazide) %>% summarise(adherence=mean(adherent)) %>% 
  ggplot(data=., aes(x=interval, y=adherence, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(1,24)) + ylim(c(0,1)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('Proportion') +
  xlab('') + theme(text = element_text(size = 10))

######################################################
######## Abnormal BP balance for total effect ########
######################################################

total_abnormalBP_1 <- sim_data[[1]] %>% group_by(interval, trt_thiazide) %>% summarise(abnormalBP=mean(abnormal_BP)) %>% 
  group_by(trt_thiazide) %>%  
  ggplot(data=., aes(x=interval, y=abnormalBP, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(1,24)) + ylim(c(0,1)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('Proportion') +
  xlab('') + theme(text = element_text(size = 10))

####################################################################################################

total_abnormalBP_2 <- sim_data[[2]] %>% group_by(interval, trt_thiazide) %>% summarise(abnormalBP=mean(abnormal_BP)) %>% 
  group_by(trt_thiazide) %>%  
  ggplot(data=., aes(x=interval, y=abnormalBP, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(1,24)) + ylim(c(0,1)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('Proportion') +
  xlab('') + theme(text = element_text(size = 10))

####################################################################################################

total_abnormalBP_3 <- sim_data[[3]] %>% group_by(interval, trt_thiazide) %>% summarise(abnormalBP=mean(abnormal_BP)) %>% 
  group_by(trt_thiazide) %>%  
  ggplot(data=., aes(x=interval, y=abnormalBP, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(1,24)) + ylim(c(0,1)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('Proportion') +
  xlab('') + theme(text = element_text(size = 10))

######################################################
############ AKI balance for total effect ############
######################################################

total_AKI_1 <- sim_data[[1]] %>% group_by(interval, trt_thiazide) %>% summarise(AKI=mean(AKI)) %>% 
  group_by(trt_thiazide) %>%  
  ggplot(data=., aes(x=interval, y=AKI, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(1,24)) + ylim(c(0,0.1)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('Proportion') +
  xlab('Month') + theme(text = element_text(size = 10))

####################################################################################################

total_AKI_2 <- sim_data[[2]] %>% group_by(interval, trt_thiazide) %>% summarise(AKI=mean(AKI)) %>% 
  group_by(trt_thiazide) %>%  
  ggplot(data=., aes(x=interval, y=AKI, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(1,24)) + ylim(c(0,0.1)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('Proportion') +
  xlab('Month') + theme(text = element_text(size = 10))

####################################################################################################

total_AKI_3 <- sim_data[[3]] %>% group_by(interval, trt_thiazide) %>% summarise(AKI=mean(AKI)) %>% 
  group_by(trt_thiazide) %>%  
  ggplot(data=., aes(x=interval, y=AKI, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(1,24)) + ylim(c(0,0.1)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('Proportion') +
  xlab('Month') + theme(text = element_text(size = 10))


######################################################
################# Separable effect  ##################
######################################################

sim_data[[1]]$Wt1_num <- predict(lm(dead ~ abnormal_BP + adherent + AKI, data=sim_data[[1]]), 
                                 newdata=sim_data[[1]] %>% mutate(trt_thiazide=0))

sim_data[[1]]$Wt1_denom <- predict(lm(dead ~ abnormal_BP + adherent + AKI, data=sim_data[[1]]), 
                                   newdata=sim_data[[1]] %>% mutate(trt_thiazide=1))

sim_data[[1]]$Wt2_num <- predict(lm(abnormal_BP ~ lag_adherent + lag_adherent:trt_thiazide, data=sim_data[[1]]), 
                                 newdata=sim_data[[1]] %>% mutate(trt_thiazide=0))

sim_data[[1]]$Wt2_denom <- predict(lm(abnormal_BP ~ lag_adherent + lag_adherent:trt_thiazide, data=sim_data[[1]]), 
                                   newdata=sim_data[[1]] %>% mutate(trt_thiazide=1))

sim_data[[1]] <- sim_data[[1]] %>% 
  mutate(Wt2_num = case_when(abnormal_BP == 1 ~ Wt2_num, abnormal_BP == 0 ~ 1-Wt2_num),
         Wt2_denom = case_when(abnormal_BP == 1 ~ Wt2_denom, abnormal_BP == 0 ~ 1-Wt2_denom)) %>%
  group_by(id) %>%
  mutate(Wt1 = cumprod(Wt1_num/Wt1_denom),
         Wt2 = cumprod(Wt2_num/Wt2_denom)) %>%
  mutate(Wt=Wt1*Wt2) %>%
  ungroup()

separable_dead_1 <- bind_rows(
  data.frame(
    sim_data[[1]] %>% filter(trt_thiazide==1) %>% group_by(interval) %>% summarise(cond_risk=mean(Wt*dead)/mean(Wt)) %>% 
      mutate(CI=1-cumprod(1-cond_risk)),
    trt_thiazide = 0
  ),
  data.frame(
    sim_data[[1]] %>% filter(trt_thiazide==1) %>% group_by(interval) %>% summarise(cond_risk=mean(dead)) %>% 
      mutate(CI=1-cumprod(1-cond_risk)),
    trt_thiazide = 1
  )
) %>% 
  {bind_rows(as.data.frame(.),data.frame(interval=c(0,0),trt_thiazide=c(0,1),cond_risk=c(0,0),CI=c(0,0)))} %>%
  ggplot(data=., aes(x=interval, y=CI, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(0,24)) + ylim(c(0,0.5)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('') +
  xlab('') + theme(text = element_text(size = 10))

####################################################################################################

sim_data[[2]]$Wt1_num <- predict(lm(dead ~ abnormal_BP + adherent + AKI, data=sim_data[[2]]), 
                                 newdata=sim_data[[2]] %>% mutate(trt_thiazide=0))

sim_data[[2]]$Wt1_denom <- predict(lm(dead ~ abnormal_BP + adherent + AKI, data=sim_data[[2]]), 
                                   newdata=sim_data[[2]] %>% mutate(trt_thiazide=1))

sim_data[[2]]$Wt2_num <- predict(lm(abnormal_BP ~ lag_adherent + lag_adherent:trt_thiazide, data=sim_data[[2]]), 
                                 newdata=sim_data[[2]] %>% mutate(trt_thiazide=0))

sim_data[[2]]$Wt2_denom <- predict(lm(abnormal_BP ~ lag_adherent + lag_adherent:trt_thiazide, data=sim_data[[2]]), 
                                   newdata=sim_data[[2]] %>% mutate(trt_thiazide=1))

sim_data[[2]] <- sim_data[[2]] %>% 
  mutate(Wt2_num = case_when(abnormal_BP == 1 ~ Wt2_num, abnormal_BP == 0 ~ 1-Wt2_num),
         Wt2_denom = case_when(abnormal_BP == 1 ~ Wt2_denom, abnormal_BP == 0 ~ 1-Wt2_denom)) %>%
  group_by(id) %>%
  mutate(Wt1 = cumprod(Wt1_num/Wt1_denom),
         Wt2 = cumprod(Wt2_num/Wt2_denom)) %>%
  mutate(Wt=Wt1*Wt2) %>%
  ungroup()

separable_dead_2 <- bind_rows(
  data.frame(
    sim_data[[2]] %>% filter(trt_thiazide==1) %>% group_by(interval) %>% summarise(cond_risk=mean(Wt*dead)/mean(Wt)) %>% 
      mutate(CI=1-cumprod(1-cond_risk)),
    trt_thiazide = 0
  ),
  data.frame(
    sim_data[[2]] %>% filter(trt_thiazide==1) %>% group_by(interval) %>% summarise(cond_risk=mean(dead)) %>% 
      mutate(CI=1-cumprod(1-cond_risk)),
    trt_thiazide = 1
  )
) %>% 
  {bind_rows(as.data.frame(.),data.frame(interval=c(0,0),trt_thiazide=c(0,1),cond_risk=c(0,0),CI=c(0,0)))} %>%
  ggplot(data=., aes(x=interval, y=CI, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(0,24)) + ylim(c(0,0.5)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('') +
  xlab('') + theme(text = element_text(size = 10))

####################################################################################################

sim_data[[3]]$Wt1_num <- predict(lm(dead ~ abnormal_BP + adherent + AKI, data=sim_data[[3]]), 
                                 newdata=sim_data[[3]] %>% mutate(trt_thiazide=0))

sim_data[[3]]$Wt1_denom <- predict(lm(dead ~ abnormal_BP + adherent + AKI, data=sim_data[[3]]), 
                                   newdata=sim_data[[3]] %>% mutate(trt_thiazide=1))

sim_data[[3]]$Wt2_num <- predict(lm(abnormal_BP ~ lag_adherent + lag_adherent:trt_thiazide, data=sim_data[[3]]), 
                                 newdata=sim_data[[3]] %>% mutate(trt_thiazide=0))

sim_data[[3]]$Wt2_denom <- predict(lm(abnormal_BP ~ lag_adherent + lag_adherent:trt_thiazide, data=sim_data[[3]]), 
                                   newdata=sim_data[[3]] %>% mutate(trt_thiazide=1))

sim_data[[3]] <- sim_data[[3]] %>% 
  mutate(Wt2_num = case_when(abnormal_BP == 1 ~ Wt2_num, abnormal_BP == 0 ~ 1-Wt2_num),
         Wt2_denom = case_when(abnormal_BP == 1 ~ Wt2_denom, abnormal_BP == 0 ~ 1-Wt2_denom)) %>%
  group_by(id) %>%
  mutate(Wt1 = cumprod(Wt1_num/Wt1_denom),
         Wt2 = cumprod(Wt2_num/Wt2_denom)) %>%
  mutate(Wt=Wt1*Wt2) %>%
  ungroup()

separable_dead_3 <- bind_rows(
  data.frame(
    sim_data[[3]] %>% filter(trt_thiazide==1) %>% group_by(interval) %>% summarise(cond_risk=mean(Wt*dead)/mean(Wt)) %>% 
      mutate(CI=1-cumprod(1-cond_risk)),
    trt_thiazide = 0
  ),
  data.frame(
    sim_data[[3]] %>% filter(trt_thiazide==1) %>% group_by(interval) %>% summarise(cond_risk=mean(dead)) %>% 
      mutate(CI=1-cumprod(1-cond_risk)),
    trt_thiazide = 1
  )
) %>% 
  {bind_rows(as.data.frame(.),data.frame(interval=c(0,0),trt_thiazide=c(0,1),cond_risk=c(0,0),CI=c(0,0)))} %>%
  ggplot(data=., aes(x=interval, y=CI, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(0,24)) + ylim(c(0,0.5)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('') +
  xlab('') + theme(text = element_text(size = 10))

######################################################
############## Check separable effect  ###############
######################################################

separable_dead_1_check <- bind_rows(
  data.frame(
    sim_data[[4]] %>% filter(trt_ZA==1 & trt_ZY==0) %>% group_by(interval) %>% summarise(cond_risk=mean(dead)) %>% 
      mutate(CI=1-cumprod(1-cond_risk)),
    trt_thiazide = 0
  ),
  data.frame(
    sim_data[[1]] %>% filter(trt_thiazide==1) %>% group_by(interval) %>% summarise(cond_risk=mean(dead)) %>% 
      mutate(CI=1-cumprod(1-cond_risk)),
    trt_thiazide = 1
  )
) %>% 
  {bind_rows(as.data.frame(.),data.frame(interval=c(0,0),trt_thiazide=c(0,1),cond_risk=c(0,0),CI=c(0,0)))} %>%
  ggplot(data=., aes(x=interval, y=CI, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(0,24)) + ylim(c(0,0.5)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('') +
  xlab('') + theme(text = element_text(size = 10))

####################################################################################################

separable_dead_2_check <- bind_rows(
  data.frame(
    sim_data[[5]] %>% filter(trt_ZA==1 & trt_ZY==0) %>% group_by(interval) %>% summarise(cond_risk=mean(dead)) %>% 
      mutate(CI=1-cumprod(1-cond_risk)),
    trt_thiazide = 0
  ),
  data.frame(
    sim_data[[2]] %>% filter(trt_thiazide==1) %>% group_by(interval) %>% summarise(cond_risk=mean(dead)) %>% 
      mutate(CI=1-cumprod(1-cond_risk)),
    trt_thiazide = 1
  )
) %>% 
  {bind_rows(as.data.frame(.),data.frame(interval=c(0,0),trt_thiazide=c(0,1),cond_risk=c(0,0),CI=c(0,0)))} %>%
  ggplot(data=., aes(x=interval, y=CI, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(0,24)) + ylim(c(0,0.5)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('') +
  xlab('') + theme(text = element_text(size = 10))

####################################################################################################

separable_dead_3_check <- bind_rows(
  data.frame(
    sim_data[[6]] %>% filter(trt_ZA==1 & trt_ZY==0) %>% group_by(interval) %>% summarise(cond_risk=mean(dead)) %>% 
      mutate(CI=1-cumprod(1-cond_risk)),
    trt_thiazide = 0
  ),
  data.frame(
    sim_data[[3]] %>% filter(trt_thiazide==1) %>% group_by(interval) %>% summarise(cond_risk=mean(dead)) %>% 
      mutate(CI=1-cumprod(1-cond_risk)),
    trt_thiazide = 1
  )
) %>% 
  {bind_rows(as.data.frame(.),data.frame(interval=c(0,0),trt_thiazide=c(0,1),cond_risk=c(0,0),CI=c(0,0)))} %>%
  ggplot(data=., aes(x=interval, y=CI, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(0,24)) + ylim(c(0,0.5)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('') +
  xlab('') + theme(text = element_text(size = 10))

######################################################
####### Adherence balance for separable effect #######
######################################################

separable_adherence_1 <- bind_rows(
  data.frame(
    sim_data[[4]] %>% filter(trt_ZA==1 & trt_ZY==0) %>% group_by(interval) %>% summarise(adherence=mean(adherent)),
    trt_thiazide = 0
  ),
  data.frame(
    sim_data[[1]] %>% filter(trt_thiazide==1) %>% group_by(interval) %>% summarise(adherence=mean(adherent)),
    trt_thiazide = 1
  )
) %>% 
  ggplot(data=., aes(x=interval, y=adherence, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(1,24)) + ylim(c(0,1)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('') +
  xlab('') + theme(text = element_text(size = 10))

####################################################################################################

separable_adherence_2 <- bind_rows(
  data.frame(
    sim_data[[5]] %>% filter(trt_ZA==1 & trt_ZY==0) %>% group_by(interval) %>% summarise(adherence=mean(adherent)),
    trt_thiazide = 0
  ),
  data.frame(
    sim_data[[2]] %>% filter(trt_thiazide==1) %>% group_by(interval) %>% summarise(adherence=mean(adherent)),
    trt_thiazide = 1
  )
) %>% 
  ggplot(data=., aes(x=interval, y=adherence, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(1,24)) + ylim(c(0,1)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('') +
  xlab('') + theme(text = element_text(size = 10))

####################################################################################################

separable_adherence_3 <- bind_rows(
  data.frame(
    sim_data[[6]] %>% filter(trt_ZA==1 & trt_ZY==0) %>% group_by(interval) %>% summarise(adherence=mean(adherent)),
    trt_thiazide = 0
  ),
  data.frame(
    sim_data[[3]] %>% filter(trt_thiazide==1) %>% group_by(interval) %>% summarise(adherence=mean(adherent)),
    trt_thiazide = 1
  )
) %>% 
  ggplot(data=., aes(x=interval, y=adherence, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(1,24)) + ylim(c(0,1)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('') +
  xlab('') + theme(text = element_text(size = 10))

######################################################
###### Abnormal BP balance for separable effect ######
######################################################

separable_abnormalBP_1 <- bind_rows(
  data.frame(
    sim_data[[4]] %>% filter(trt_ZA==1 & trt_ZY==0) %>% group_by(interval) %>% summarise(abnormalBP=mean(abnormal_BP)),
    trt_thiazide = 0
  ),
  data.frame(
    sim_data[[1]] %>% filter(trt_thiazide==1) %>% group_by(interval) %>% summarise(abnormalBP=mean(abnormal_BP)),
    trt_thiazide = 1
  )
) %>% 
  ggplot(data=., aes(x=interval, y=abnormalBP, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(1,24)) + ylim(c(0,1)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('') +
  xlab('') + theme(text = element_text(size = 10))

####################################################################################################

separable_abnormalBP_2 <- bind_rows(
  data.frame(
    sim_data[[5]] %>% filter(trt_ZA==1 & trt_ZY==0) %>% group_by(interval) %>% summarise(abnormalBP=mean(abnormal_BP)),
    trt_thiazide = 0
  ),
  data.frame(
    sim_data[[2]] %>% filter(trt_thiazide==1) %>% group_by(interval) %>% summarise(abnormalBP=mean(abnormal_BP)),
    trt_thiazide = 1
  )
) %>% 
  ggplot(data=., aes(x=interval, y=abnormalBP, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(1,24)) + ylim(c(0,1)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('') +
  xlab('') + theme(text = element_text(size = 10))

####################################################################################################

separable_abnormalBP_3 <- bind_rows(
  data.frame(
    sim_data[[6]] %>% filter(trt_ZA==1 & trt_ZY==0) %>% group_by(interval) %>% summarise(abnormalBP=mean(abnormal_BP)),
    trt_thiazide = 0
  ),
  data.frame(
    sim_data[[3]] %>% filter(trt_thiazide==1) %>% group_by(interval) %>% summarise(abnormalBP=mean(abnormal_BP)),
    trt_thiazide = 1
  )
) %>% 
  ggplot(data=., aes(x=interval, y=abnormalBP, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(1,24)) + ylim(c(0,1)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('') +
  xlab('') + theme(text = element_text(size = 10))

######################################################
########## AKI balance for separable effect ##########
######################################################

separable_AKI_1 <- bind_rows(
  data.frame(
    sim_data[[4]] %>% filter(trt_ZA==1 & trt_ZY==0) %>% group_by(interval) %>% summarise(AKI=mean(AKI)),
    trt_thiazide = 0
  ),
  data.frame(
    sim_data[[1]] %>% filter(trt_thiazide==1) %>% group_by(interval) %>% summarise(AKI=mean(AKI)),
    trt_thiazide = 1
  )
) %>% 
  ggplot(data=., aes(x=interval, y=AKI, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(1,24)) + ylim(c(0,0.1)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('') +
  xlab('Month') + theme(text = element_text(size = 10))

####################################################################################################

separable_AKI_2 <- bind_rows(
  data.frame(
    sim_data[[5]] %>% filter(trt_ZA==1 & trt_ZY==0) %>% group_by(interval) %>% summarise(AKI=mean(AKI)),
    trt_thiazide = 0
  ),
  data.frame(
    sim_data[[2]] %>% filter(trt_thiazide==1) %>% group_by(interval) %>% summarise(AKI=mean(AKI)),
    trt_thiazide = 1
  )
) %>% 
  ggplot(data=., aes(x=interval, y=AKI, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(1,24)) + ylim(c(0,0.1)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('') +
  xlab('Month') + theme(text = element_text(size = 10))

####################################################################################################

separable_AKI_3 <- bind_rows(
  data.frame(
    sim_data[[6]] %>% filter(trt_ZA==1 & trt_ZY==0) %>% group_by(interval) %>% summarise(AKI=mean(AKI)),
    trt_thiazide = 0
  ),
  data.frame(
    sim_data[[3]] %>% filter(trt_thiazide==1) %>% group_by(interval) %>% summarise(AKI=mean(AKI)),
    trt_thiazide = 1
  )
) %>% 
  ggplot(data=., aes(x=interval, y=AKI, col=as.factor(trt_thiazide), linetype=as.factor(trt_thiazide))) +
  geom_point(size=0.25) + geom_line(size=0.25) + xlim(c(1,24)) + ylim(c(0,0.1)) + theme_classic() +
  theme(legend.title=element_blank()) +
  scale_color_hue(labels=c('ACEI','Thiazide')) +
  scale_linetype(labels=c('ACEI','Thiazide')) +
  ylab('') +
  xlab('Month') + theme(text = element_text(size = 10))

######################################################
######################## Graph #######################
######################################################

####################################################################################################

adherence_model_1_1 <- ggarrange(total_dead_1, separable_dead_1,
                                 ncol=2, nrow=1,
                                 common.legend=T,
                                 legend='none') %>% 
  annotate_figure(., top = text_grob("Cumulative incidence of death", color = "black", size = 12))

adherence_model_1_2 <- ggarrange(total_adherence_1, separable_adherence_1,
                                 ncol=2, nrow=1,
                                 common.legend=T,
                                 legend='none') %>%
  annotate_figure(., top = text_grob("Probability of adherence", color = "black", size = 12))

adherence_model_1_3 <- ggarrange(total_abnormalBP_1, separable_abnormalBP_1,
                                 ncol=2, nrow=1,
                                 common.legend=T,
                                 legend='none') %>%
  annotate_figure(., top = text_grob("Probability of abnormal BP", color = "black", size = 12))

adherence_model_1_4 <- ggarrange(total_AKI_1, separable_AKI_1,
                                 ncol=2, nrow=1,
                                 common.legend=T,
                                 legend='bottom') %>%
  annotate_figure(., top = text_grob("Probability of AKI", color = "black", size = 12))

ggarrange(adherence_model_1_1,
          adherence_model_1_2,
          adherence_model_1_3,
          adherence_model_1_4, 
          ncol=1, nrow=4) %>%
  annotate_figure(., top = text_grob("Adherence model 1", color = "black", 
                                     face='bold', hjust=1.65, size = 18),
                  bottom = text_grob("Total effect                                                                Separable effect", 
                                     color = "Black", 
                                     hjust=0.43, size = 12))

ggsave("./figures/adherence_model_1.png", 
       dpi=600, bg = "white", height = 225, width = 200, units = "mm")

####################################################################################################

adherence_model_2_1 <- ggarrange(total_dead_2, separable_dead_2,
                                 ncol=2, nrow=1,
                                 common.legend=T,
                                 legend='none') %>% 
  annotate_figure(., top = text_grob("Cumulative incidence of death", color = "black", size = 12))

adherence_model_2_2 <- ggarrange(total_adherence_2, separable_adherence_2,
                                 ncol=2, nrow=1,
                                 common.legend=T,
                                 legend='none') %>%
  annotate_figure(., top = text_grob("Probability of adherence", color = "black", size = 12))

adherence_model_2_3 <- ggarrange(total_abnormalBP_2, separable_abnormalBP_2,
                                 ncol=2, nrow=1,
                                 common.legend=T,
                                 legend='none') %>%
  annotate_figure(., top = text_grob("Probability of abnormal BP", color = "black", size = 12))

adherence_model_2_4 <- ggarrange(total_AKI_2, separable_AKI_2,
                                 ncol=2, nrow=1,
                                 common.legend=T,
                                 legend='bottom') %>%
  annotate_figure(., top = text_grob("Probability of AKI", color = "black", size = 12))

ggarrange(adherence_model_2_1,
          adherence_model_2_2,
          adherence_model_2_3,
          adherence_model_2_4, 
          ncol=1, nrow=4) %>%
  annotate_figure(., top = text_grob("Adherence model 2", color = "black", 
                                     face='bold', hjust=1.65, size = 18),
                  bottom = text_grob("Total effect                                                                Separable effect", 
                                     color = "Black", 
                                     hjust=0.43, size = 12))

ggsave("./figures/adherence_model_2.png", 
       dpi=600, bg = "white", height = 225, width = 200, units = "mm")

####################################################################################################

adherence_model_3_1 <- ggarrange(total_dead_3, separable_dead_3,
                                 ncol=2, nrow=1,
                                 common.legend=T,
                                 legend='none') %>% 
  annotate_figure(., top = text_grob("Cumulative incidence of death", color = "black", size = 12))

adherence_model_3_2 <- ggarrange(total_adherence_3, separable_adherence_3,
                                 ncol=2, nrow=1,
                                 common.legend=T,
                                 legend='none') %>%
  annotate_figure(., top = text_grob("Probability of adherence", color = "black", size = 12))

adherence_model_3_3 <- ggarrange(total_abnormalBP_3, separable_abnormalBP_3,
                                 ncol=2, nrow=1,
                                 common.legend=T,
                                 legend='none') %>%
  annotate_figure(., top = text_grob("Probability of abnormal BP", color = "black", size = 12))

adherence_model_3_4 <- ggarrange(total_AKI_3, separable_AKI_3,
                                 ncol=2, nrow=1,
                                 common.legend=T,
                                 legend='bottom') %>%
  annotate_figure(., top = text_grob("Probability of AKI", color = "black", size = 12))

ggarrange(adherence_model_3_1,
          adherence_model_3_2,
          adherence_model_3_3,
          adherence_model_3_4, 
          ncol=1, nrow=4) %>%
  annotate_figure(., top = text_grob("Adherence model 3", color = "black", 
                                     face='bold', hjust=1.65, size = 18),
                  bottom = text_grob("Total effect                                                                Separable effect", 
                                     color = "Black", 
                                     hjust=0.43, size = 12))

ggsave("./figures/adherence_model_3.png", 
       dpi=600, bg = "white", height = 225, width = 200, units = "mm")