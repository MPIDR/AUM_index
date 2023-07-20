##--------------------------------------------------##
##
## Computing the AUM from observed period death rates
##
## Generating Figures 1, 2 and 3 of the paper
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
  save(DF,file="out/F1-F3.Rdata")
  
}else{
  load("out/F1-F3.Rdata")
}


##---- Figure 1 -----
DFsub <- DF %>% filter(Age==0,(Sex=="male" & country=="ITA")|
                         (Sex=="female" & country=="SWE"))

## left panel - lifespan variability -- color ex
f1A <- DFsub %>%
  select(Year,Sex,ex,edag,sd,country) %>% 
  pivot_longer(cols=c(edag,sd)) %>% 
  mutate(name2 = paste(Sex,name,sep="-"),
         Sex2 =  recode(Sex, female = "Swedish females" , male = "Italian males" ),
         Sex3 = factor(Sex2,levels = c("Swedish females","Italian males"))) %>% 
  ggplot(aes(x=Year,y=value,color=ex,shape=name))+
  geom_point(size=2.15) +
  facet_wrap(Sex3~.,ncol = 1) +
  annotate("text", x=1975, y=25, label= expression(sigma[T]),
           col="black",size=7) +
  annotate("text", x=1905, y=15, label= expression(e),
           col="black",size=7) +
  scale_shape(guide = 'none') +
  theme_bw(base_size = 16) +
  labs(x="Year",y="Lifespan variability") +
  scale_color_viridis(name=expression(e[0]),option = "A") 
f1A

## right panel - AUM
## color ex
f1B <- DFsub %>%
  select(Year,Sex,ex,AUM,country) %>% 
  mutate(Sex2 =  recode(Sex, female = "Swedish females" , male = "Italian males" ),
         Sex3 = factor(Sex2,levels = c("Swedish females","Italian males"))) %>% 
  ggplot(aes(x=Year,y=AUM,color=ex))+
  geom_point(size=2.15) +
  facet_wrap(Sex3~.,ncol = 1) +
  scale_color_viridis(name=expression(e[0]),option = "A") +
  theme_bw(base_size = 16) +
  labs(x="Year",y="AUM")
f1B

f1A + f1B + plot_layout(guides = 'collect')

ggsave("out/F1-pre.pdf",width = 12,height = 8)

##---- Figure 2 -----
DF %>% filter(Age==0) %>% 
  ggplot(aes(x=Year,y=AUM,color=Sex))+
  geom_point(size=1.05) +
  annotate("text", x=1890, y=0.665, label= "females",
           col="#998EC3",size=10,fontface="italic") +
  annotate("text", x=1985, y=0.86, label= "males",
           col="#F1A340",size=10,fontface="italic") +
  scale_color_manual(values=c("#998EC3","#F1A340")) +
  theme_bw(base_size = 16) +
  labs(x="Year",y="AUM") +
  theme(legend.position = "none")


ggsave("out/F2.pdf",width = 12,height = 8)

##---- Figure 2b -----
f2bA <- DF %>% filter(Age==0, Sex=="female") %>% 
  ggplot(aes(x=edag,y=AUM,color=Year))+
  geom_point() +
  theme_bw(base_size = 16) +
  labs(x="e",y="AUM") +
  scale_color_viridis(option = "C")
f2bA

f2bB <- DF %>% filter(Age==0, Sex=="female") %>% 
  mutate(H=edag/ex) %>% 
  ggplot(aes(x=H,y=AUM,color=Year))+
  geom_point() +
  theme_bw(base_size = 16) +
  labs(x="Entropy",y=NULL) +
  scale_color_viridis(option = "C")
f2bB

f2bA + f2bB + plot_layout(guides = 'collect')

ggsave("out/F2b-pre.pdf",width = 12,height = 8)


##---- Figure 3 -----
DF2 <- DF %>% 
  mutate(Sex=recode(Sex,male="Males",female="Females"))

DF2 %>% filter(Age==0) %>% 
  ggplot(aes(x=ex,y=AUM,color=Year))+
  geom_point()+geom_smooth(col="black", fill = NA)+
  facet_wrap(.~Sex)+
  theme_bw(base_size = 16) +
  labs(x="Period Life Expectancy",y="AUM") +
  scale_color_viridis(option = "C") 
ggsave("out/F3.pdf",width = 12,height = 8)





