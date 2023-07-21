##--------------------------------------------------##
##
## Computing the AUM from observed period death rates
## and modelled period life-table death rates
##
## Generating Figure 5 and Table 1 of the paper
##
##--------------------------------------------------##

## cleaning the workspace
rm(list=ls(all=TRUE))

## loading packages
library(HMDHFDplus)
library(tidyverse)
library(patchwork)
library(viridis)

## loading previous results?
UPDATE <- FALSE

if (UPDATE){

  ## loading functions
  source("fun/funs-AUM.R")
  source("fun/funs-ineq.R")
  
  ## insert your HMD username and password
  myusername <- myusername 
  mypassword <- mypassword
  
  ## select country
  country <- c( "SWE")
  cat("analysing",country,cat="\n")
  
  ## observed death rates
  DEA <- readHMDweb(CNTRY=country,item = "Deaths_1x1",
                    username=myusername, password=mypassword, fixup=T)
  DEAf <- DEA %>% 
    select(Year,Age,Deaths=Female)
  
  EXPO <- readHMDweb(CNTRY=country,item = "Exposures_1x1",
                     username=myusername, password=mypassword, fixup=T)
  
  EXPOf <- EXPO %>% 
    select(Year,Age,Exposures=Female)
  
  RATES <- DEAf %>% 
    left_join(EXPOf,by = join_by(Year, Age)) %>% 
    mutate(Mx=Deaths/Exposures)
  
  ## life table death rates
  LTf <- readHMDweb(CNTRY=country,item = "fltper_1x1",
                    username=myusername, password=mypassword, fixup=T)
  
  ## compute AUM for observed death rates
  DF1 <- RATES %>% group_by(Year) %>% 
    mutate(edag=CovTLambdaTCond(age=Age,mx=Mx),
           sd=sqrt(varTCond(age=Age,mx=Mx)),
           AUM=edag/sd,
           type = "observed rates")
  
  ## compute AUM for life-table death rates
  DF2 <- LTf %>% group_by(Year) %>% 
    mutate(edag=CovTLambdaTCond(age=Age,mx=mx),
           sd=sqrt(varTCond(age=Age,mx=mx)),
           AUM=edag/sd,
           type = "life-table rates")
  
  ## compute AUM using life-table functions
  DF3 <- LTf %>% group_by(Year) %>%
    mutate(edag=edag_fun(age=Age,mx=mx,dx=dx,lx=lx,ex=ex,ax=ax),
           sd=sd_fun(age=Age,mx=mx,dx=dx,lx=lx,ex=ex,ax=ax),
           AUM=edag/sd,
           type = "life-table functions")
  
  ## combine
  DF <- DF1 %>% bind_rows(DF2) %>% bind_rows(DF3)
  
  ## save results
  save(DF,file="out/F6.Rdata")
  
}else{
  load("out/F6.Rdata")
}


##---- Figure 1 Supp Material -----
DF %>% 
  filter(type!="observed rates") %>%
    mutate(type=factor(type,levels = c("life-table rates","life-table functions")),
           type=recode(type, "life-table rates"="own routines",
                       "life-table functions"="conventional routines")) %>% 
    ggplot(aes(x=Age,y=AUM,color=Year))+
    geom_point() +
    facet_wrap(.~type) +
    scale_color_viridis() +
    theme_bw(base_size = 16) 

ggsave("out/SM-F1.pdf",width = 12,height = 8)


##---- Table 1 -----
## computing the percentage relative difference
DF <- DF2 %>% bind_rows(DF3) %>%
  select(Year,Age,edag,sd,AUM,type) %>%
  mutate(type=recode(type, 'life-table rates' = "true" , 'life-table functions' = "bias" )) %>% 
  pivot_wider(names_from = "type",values_from = c(edag,sd,AUM)) %>% 
  mutate(bias_AUM=100*(AUM_bias-AUM_true)/AUM_true,
         bias_edag=100*(edag_bias-edag_true)/edag_true,
         bias_sd=100*(sd_bias-sd_true)/sd_true)

DF %>% 
  ggplot(aes(x=Age,y=bias_sd,color=Year))+
  geom_point() +
  scale_color_viridis() +
  theme_bw(base_size = 16) 

## summarizing by years and ages
DFsumm <- DF %>% ungroup() %>% 
  group_by(Age) %>% 
  summarise(AUMe = mean(bias_edag,na.rm=T),
            AUMs = mean(bias_sd,na.rm=T),
            AUMb = mean(bias_AUM,na.rm=T)) %>% 
  mutate(AgeGroup = case_when(
    Age == 0 ~ 1,
    Age > 0 & Age <=60 ~ 2,
    Age > 60 & Age <=80 ~ 3,
    Age > 80 & Age <=90 ~ 4,
    Age > 90 & Age <=100 ~ 5,
    Age > 100 & Age <=105 ~ 6,
    Age > 105 & Age <=110 ~ 7)) %>% 
  ungroup() %>% group_by(AgeGroup) %>% 
  summarise(AUMe = mean(AUMe,na.rm=T),
            AUMs = mean(AUMs,na.rm=T),
            AUMb = mean(AUMb,na.rm=T))

## table 1
print(round(DFsumm,2),sep="&")



