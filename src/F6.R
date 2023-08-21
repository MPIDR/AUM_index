##--------------------------------------------------##
##
## Computing the AUM from observed period death rates
## and modelled period life-table death rates
##
## Generating Figure 6 of the paper
##
##  sessionInfo() details:
##
## R version 4.2.3 (2023-03-15 ucrt)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 19044)
##
## attached base packages:
## stats     graphics  grDevices utils     datasets  methods   base     
##
## other attached packages:
## viridis_0.6.3     viridisLite_0.4.2 patchwork_1.1.2   lubridate_1.9.2  
## forcats_1.0.0     stringr_1.5.0     dplyr_1.1.2       purrr_1.0.1      
## readr_2.1.4       tidyr_1.3.0       tibble_3.2.1      ggplot2_3.4.2    
## tidyverse_2.0.0   HMDHFDplus_2.0.1 
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


##---- Figure 6 -----
DF %>% 
  filter(type!="life-table functions") %>% 
    mutate(type=factor(type,levels = c("observed rates","life-table rates"))) %>% 
    ggplot(aes(x=Age,y=AUM,color=Year))+
    geom_point() +
    facet_wrap(.~type) +
    scale_color_viridis() +
    theme_bw(base_size = 16) 

ggsave("out/F6.pdf",width = 12,height = 8)


