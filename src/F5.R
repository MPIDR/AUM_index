##--------------------------------------------------##
##
## Computing the AUM from observed cohort death rates
##
## Generating Figure 5 of the paper
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
  
  ## insert your HMD username and password
  myusername <- myusername 
  mypassword <- mypassword
  
  ## Get USA deaths 1x1 (by single year of age and year)
  countries <- c( "FRATNP", "ITA", "SWE", "FIN", "DNK","NLD")
  n_cou <- length(countries)
  
  jj <- 1
  for (jj in 1:n_cou){
    cat("analysing",countries[jj]," - country",jj,"/",n_cou,cat="\n")
    LTf <- readHMDweb(CNTRY=countries[jj],item = "fltcoh_1x1",
                      username=myusername, password=mypassword, fixup=T)
    LTm <- readHMDweb(CNTRY=countries[jj],item = "mltcoh_1x1",
                      username=myusername, password=mypassword, fixup=T)
    
    LTf <- LTf %>% mutate(sex="female")
    LTm <- LTm %>% mutate(sex="male")
    LT <- LTf %>% bind_rows(LTm)
    
    ## ages and years
    ages <- unique(LT$Age)
    years <- unique(LT$Year)
    m <- length(ages)
    n <- length(years)
    
    ## compute AUM
    DF.temp <- LT %>% group_by(Year,sex) %>% 
      mutate(ex=ETCond(age=Age,mx=mx),
             edag=CovTLambdaTCond(age=Age,mx=mx),
             sd=sqrt(varTCond(age=Age,mx=mx)),
             AUM=edag/sd,
             country=countries[jj])
    
    ## 
    if (jj==1){
      DF <- DF.temp
    }else{
      DF <- DF %>% bind_rows(DF.temp)
    }
    
  }
  
  ## save results
  save(DF,file="out/F5.Rdata")

}else{
  
  load("out/F5.Rdata")
  
}


##---- Figure 5 -----
DF2 <- DF %>% 
  mutate(sex=recode(sex,male="Males",female="Females"))

DF2 %>% filter(Age==0) %>% 
  ggplot(aes(x=ex,y=AUM,color=Year))+
  geom_point()+geom_smooth(col="black", fill = NA)+
  facet_wrap(.~sex)+
  theme_bw(base_size = 16) +
  labs(x="Cohort Life Expectancy",y="AUM") +
  scale_color_viridis(name="Cohort",option = "C") 

ggsave("out/F5.pdf",width = 12,height = 8)



