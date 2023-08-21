##--------------------------------------------------##
##
## Computing the AUM from observed cohort death rates
##
## Generating Figure 5 of the paper
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



