##--------------------------------------------------##
##
## Computing the AUM from observed period death rates
##
## Generating Figures 1 and 2 of the Supp Materials
##
##  sessionInfo() details:
##
## R version 4.3.2 (2023-10-31 ucrt)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 19045)
##
## attached base packages:
## stats     graphics  grDevices utils     datasets  methods   base     
##
## other attached packages:
## viridis_0.6.5     viridisLite_0.4.2 patchwork_1.2.0   lubridate_1.9.3  
## forcats_1.0.0     stringr_1.5.1     dplyr_1.1.4       purrr_1.0.2      
## readr_2.1.5       tidyr_1.3.1       tibble_3.2.1      ggplot2_3.4.4    
## tidyverse_2.0.0   HMDHFDplus_2.0.3 
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
  
  ## all HMD countries
  countries <- c( "AUS", "AUT", "BEL" , "BGR" ,    "BLR"   ,   "CAN",
                  "CHL", "HRV",  "HKG"    ,  "CHE"  ,     "CZE"  ,   
                  "DEUTNP", "DNK"     , "ESP"   ,    "EST"  ,   "FIN",
                  "FRATNP",  "GRC",      "HUN"     ,   "IRL" ,    "ISL",
                  "ISR", "ITA", "JPN"     , "KOR"   ,    "LTU"   ,   "LUX",
                  "LVA", "NLD",  "NOR"    , "NZL_NP",
                  "POL", "PRT", "RUS"     , "SVK"   ,    "SVN"   ,   "SWE",
                  "TWN", "UKR", "GBR_NP" , 
                  "USA")
  n_cou <- length(countries)
  
  jj <- 1
  for (jj in 1:n_cou){
    cat("analysing",countries[jj]," - country",jj,"/",n_cou,cat="\n")
    DEA <- readHMDweb(CNTRY=countries[jj],item = "Deaths_1x1",
                      username=myusername, password=mypassword, fixup=T)
    DEAf <- DEA %>% 
      select(Year,Age,Deaths=Female) %>% 
      mutate(Sex="female")
    DEAm <- DEA %>% 
      select(Year,Age,Deaths=Male) %>% 
      mutate(Sex="male")
    DEAall <- DEAf %>% 
      bind_rows(DEAm)
    
    EXPO <- readHMDweb(CNTRY=countries[jj],item = "Exposures_1x1",
                       username=myusername, password=mypassword, fixup=T)
    
    EXPOf <- EXPO %>% 
      select(Year,Age,Exposures=Female) %>% 
      mutate(Sex="female")
    EXPOm <- EXPO %>% 
      select(Year,Age,Exposures=Male) %>% 
      mutate(Sex="male")
    EXPOall <- EXPOf %>% 
      bind_rows(EXPOm)
    
    RATES <- DEAall %>% 
      left_join(EXPOall,by = join_by(Year, Age, Sex)) %>% 
      mutate(Mx=Deaths/Exposures)
    
    ## remove 1914- 1918 in Belgium
    if (countries[jj] == "BEL"){
      RATES <- RATES %>% 
        filter(!(Year %in% c(1914:1918)))
    }
    
    ## compute AUM
    DF.temp <- RATES %>% group_by(Year,Sex) %>% 
      mutate(ex=ETCond(age=Age,mx=Mx),
             edag=CovTLambdaTCond(age=Age,mx=Mx),
             sd=sqrt(varTCond(age=Age,mx=Mx)),
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
  save(DF,file="out/F1-F4.Rdata")
  
}else{
  load("out/F1-F4.Rdata")
}


##---- Figure S1 -----
countries <- c( "FRATNP", "ITA", "SWE", "FIN", "DNK","NLD",
                "ISL", "NOR", "CHE", "GBRTENW")
DF2 <- DF %>% 
  filter(country%in%countries) %>% 
  mutate(Sex=recode(Sex,male="Males",female="Females"))

DF2 %>% filter(Age==0) %>% 
  ggplot(aes(x=ex,y=AUM,color=Year))+
  geom_point()+geom_smooth(col="black", fill = NA)+
  facet_wrap(.~Sex)+
  theme_bw(base_size = 16) +
  labs(x="Period Life Expectancy",y="AUM") +
  scale_color_viridis(option = "C") 
ggsave("out/SM-F1.pdf",width = 12,height = 8)


##---- Figure S2 -----

DF %>% filter(Sex=="female") %>% 
  filter(Age %in% seq(0,90,10)) %>% 
  mutate(ex=ex-Age) %>% 
  ## remove ISL in 1852 and 1878 (small numbers creating issues)
  filter(!(country=="ISL"&Year%in%c(1852,1878))) %>% 
  ggplot(aes(x=ex,y=AUM,color=Year))+
  geom_point()+geom_smooth(col="black", fill = NA)+
  facet_wrap(.~Age)+
  theme_bw(base_size = 16) +
  labs(x="Remaining Life Expectancy",y="AUM") +
  scale_color_viridis(option = "C") 
ggsave("out/SM-F2.pdf",width = 12,height = 8)



