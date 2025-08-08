
library(reshape2)
library(tidyverse)
library(rjags)
library(R2jags)
library(runjags)
library(tidyverse)
######Calculate turnover using abundance########

#####Monthly Abundance Data#####

allbugs.mo4<-read.csv("53bc_inputjagsdataformatted_abundance.csv")[-1]


nbugs<-length(unique(allbugs.mo4$taxon))
nyears<-length(unique(allbugs.mo4$sample.yr))
no.months<-length(unique(allbugs.mo4$sample.month))
allbugsarray.mo<-array(NA, dim=c(61,12,14))

str(allbugsarray.mo)

allbugsarray.mo[,,]<-as.array(allbugs.mo4$value)
str(allbugsarray.mo) 
bugnames<-matrix(NA, 1,61 )
str(bugnames)
bugnames[,]<-as.matrix(unique(allbugs.mo4$taxon))

###t### check a couple of values 
allbugsarray.mo[9,1,1] #should be 75; corduleg, a.9, 1993
allbugsarray.mo[40,5,1] # should be 12.5; parapsyche, e.1, 1993






#### estimate community turnover as in Shimadzu et al.
# Specify model in BUGS language
sink("coweetabugs_53_turnover.monthly.abundance.jags")
cat("
     model { 

    for(i in 1:nspp){ #  species, 61
      for(j in 1:nmonths){# months, 12
        for (k in 1:nyears){ #years,14
          N[i,j,k] ~ dnorm(nmu[i,j,k],prec)  #count, spp i, year k; each month treated as a rep for each species
          log(nmu[i,j,k]) <- logabun[i,k]  #log abundance spp i, year k
          
      }}}
        
   for(i in 1:nspp){
    for (k in 1:nyears){
   logabun[i,k] ~ dunif(-7,12)  #prior on log mean abundance
   expabun[i,k]<-exp(logabun[i,k])  ## estimated abundance
   
    }}
   
   for (k in 1:nyears){
   ttl_expected[k]<-sum(expabun[,k])  ### total for all species together, each year
   }
   
   ## calculate D1, community composition turnover
   for(i in 1:nspp){
   for (k in 1:nyears){
   p.abun[i,k]<-expabun[i,k]/ttl_expected[k]
   }}
   for(i in 1:nspp){
   for (k in 2:nyears){
   d1[i,k]<-log(p.abun[i,1]/p.abun[i,k])*p.abun[i,1]  ### first year is reference
   }}
   for (k in 2:nyears){
   D2[k]<-log(ttl_expected[k]/ttl_expected[1]) ### first year is reference
   D1[k]<--sum(d1[,k])
   D[k]<-D1[k]+D2[k]
   }
   prec <- 1/(sigma^2)
      sigma ~ dunif(0,6) 
  }  
  
    ",fill = TRUE)
sink()

# Bundle data
win.data <- list(N = allbugsarray.mo, nspp=dim(allbugsarray.mo)[1], nmonths=dim(allbugsarray.mo)[2], 
                 nyears=dim(allbugsarray.mo)[3])

# Initial values
#inits<-function(){ list(logabun.pre = rep(1,18), logabun.post =rep(1,18))}
# Initial values
ntaxa=61
nsite=1
nmonth=12
nyear=14
zst<-array(NA,dim=c(61,12,14))
for (m in 1:ntaxa){
  for(i in 1:nmonth){
    for (k in 1:nyear){
      if(is.na(allbugsarray.mo[m,i,k]))
      {zst[m,i,k]<-NA
      next}
      zst[m,i,k]<-max(allbugsarray.mo[m,i,k],na.rm=T)
    }
  }
}

inits <- function(){ list(z = zst)}

# Parameters monitored

params <- c("D2", "D1", "D","expabun")

# MCMC settings
ni <- 12000
nt <- 4
nb <- 1200
nc <- 3

# Call JAGS from R (BRT 3 min)

out.bugs.ref1<- jags(win.data, inits, params, "coweetabugs_53_turnover.monthly.abundance.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
# Summarize posteriors
print(out.bugs.ref1, dig = 2)

out.bugs.ref2<- run.jags(model= "jags_turnover_model_monhtlybugs.R",
                         monitor=params,
                         data=win.data,
                         inits=inits,
                         n.chains = nc, 
                         thin = nt, 
                         adapt=2e3,
                         sample = ni, 
                         burnin = nb,
                         summarise=F,
                         keep.jags.files = FALSE,
                        method = 'parallel'
                        )
out.bugs.summ<-summary(out.bugs.ref2)

#### these exp abundances look realistic
ref.D<-out.bugs.ref1$BUGSoutput$mean$D
ref.D.l95<-out.bugs.ref1$BUGSoutput$Lower95$D
ref.D.u95<-out.bugs.ref1$BUGSoutput$Upper95$D
ref.D1<-out.bugs.ref1$BUGSoutput$mean$D1
ref.D1.l95<-out.bugs.ref1$BUGSoutput$Lower95$D1
ref.D1.u95<-out.bugs.ref1$BUGSoutput$Upper95$D1
ref.D2<-out.bugs.ref1$BUGSoutput$mean$D2
ref.D2.l95<-out.bugs.ref1$BUGSoutput$Lower95$D2
ref.D2.u95<-out.bugs.ref1$BUGSoutput$Upper95$D2



bugstable <- matrix(rep(NA,42*4), nrow=42)#14 years x 3 D vars = 42 rows

bugstable[1,2:4]<-0#populate the first year (year 0) with 0 for no turnover
bugstable[15,2:4]<-0
bugstable[29,2:4]<-0

bugstable[2:14,2:4]<-out.bugs.summ[1:13,c(4,1,3)]#populate table with output for D
bugstable[16:28,2:4]<-out.bugs.summ[14:26,c(4,1,3)]#populate table with output for D1
bugstable[30:42,2:4]<-out.bugs.summ[27:39,c(4,1,3)]#populate table with output for D2
bugstable<-data.frame(bugstable)
bugstable[1:14,1]<-"D2"
bugstable[15:28,1]<-"D1"
bugstable[29:42,1]<-"D"
names(bugstable)[1]<-"Dvar"
names(bugstable)[2]<-"Mean"
names(bugstable)[3]<-"L95"
names(bugstable)[4]<-"U95"




year<-c(seq(1993,2006,1))
bugstable<-cbind(bugstable, year)


D.data<-cbind(year, c(0,ref.D))##0 bc jags doesn't include D0 (initial value)
D.data<-cbind(D.data, c(0, ref.D1))
D.data<-cbind(D.data, c(0, ref.D2))
D.data<-as.data.frame(D.data)
D.data<-cbind(D.data, c(0,ref.D.l95))
D.data<-cbind(D.data, c(0,ref.D.u95))
D.data<-cbind(D.data, c(0,ref.D1.l95))
D.data<-cbind(D.data, c(0,ref.D1.u95))
D.data<-cbind(D.data, c(0,ref.D2.l95))
D.data<-cbind(D.data, c(0,ref.D2.u95))
names(D.data)[2]<-"D"
names(D.data)[3]<-"D1"
names(D.data)[4]<-"D2"
names(D.data)[5]<-"D_L95"
names(D.data)[6]<-"D_U95"
names(D.data)[7]<-"D1_L95"
names(D.data)[8]<-"D1_U95"
names(D.data)[9]<-"D2_L95"
names(D.data)[10]<-"D2_U95"



#### pull out 'expected abundances'
exp.abun<-out.bugs.ref1$BUGSoutput$mean$expabun
exp.abun<-as.data.frame(exp.abun)
exp.abun.t<-t(exp.abun)  ## expected abundances, years x spp

exp.abundance<-exp.abun.t
colnames(exp.abundance)<-bugnames[1,]
#write.csv(exp.abundance, file="exp_abund_53bc_annualmeans.csv")

exp.abundance<-data.frame(exp.abundance)
exp.abundance<-cbind(exp.abundance, year)
exp.abundance$habitat<-"BC"
exp.abundance<-exp.abundance[,c(62,63,1:61)]



#write.csv(exp.abundance, file="exp_abund_53bc_abundance_annualmeans.csv")



jm.mcmc=as.mcmc(out.bugs.ref1)

nvalues<-1
year2<-c(seq(1994,2006,1))
#jm1_mcmc_combi<-as.mcmc(rbind(jm1.mcmc[[1]],jm1.mcmc[[2]], jm1.mcmc[[3]], jm1.mcmc[[4]], jm1.mcmc[[5]]  ))


plot(year2, jm.sum[65:77,"Mean"])


D2<-D.data %>%
  mutate(delD = D - lag(D, default = 0), delD1=D1 - lag(D1, default = 0),delD2=D2 - lag(D2, default = 0) )



######Rockface######
allbugs.mo4.rf<-read.csv( "53rf_inputjagsdataformatted_abundance.csv")[-1]

nbugs.rf<-length(unique(allbugs.mo4.rf$taxon))
nyears.rf<-length(unique(allbugs.mo4.rf$sample.yr))
no.months.rf<-length(unique(allbugs.mo4.rf$sample.month))


allbugsarray.mo.rf<-array(NA, dim=c(53,12,14 ))

str(allbugsarray.mo.rf)

allbugsarray.mo.rf[,,]<-as.array(allbugs.mo4.rf$value)
str(allbugsarray.mo.rf) 
bugnames.rf<-matrix(NA, 1,53 )
str(bugnames.rf)
bugnames.rf[,]<-as.matrix(unique(allbugs.mo4.rf$taxon))

###ta-da### check a couple of values 
allbugsarray.mo.rf[10,1,1] #should be 236.33, diplsp

# Specify model in BUGS language
sink("coweetabugs_53_turnover.monthly.jags")
cat("
     model { 

    for(i in 1:nspp){ #  species, 67
      for(j in 1:nmonths){# months, 12
        for (k in 1:nyears){ #years,14
          N[i,j,k] ~ dnorm(nmu[i,j,k],prec)  #count, spp i, year k; each month treated as a rep for each species
          log(nmu[i,j,k]) <- logabun[i,k]  #log abundance spp i, year k
          
      }}}
        
   for(i in 1:nspp){
    for (k in 1:nyears){
   logabun[i,k] ~ dunif(-7,12)  #prior on log mean abundance
   expabun[i,k]<-exp(logabun[i,k])  ## estimated abundance
   
    }}
   
   for (k in 1:nyears){
   ttl_expected[k]<-sum(expabun[,k])  ### total for all species together, each year
   }
   
   ## calculate D1, community composition turnover
   for(i in 1:nspp){
   for (k in 1:nyears){
   p.abun[i,k]<-expabun[i,k]/ttl_expected[k]
   }}
   for(i in 1:nspp){
   for (k in 2:nyears){
   d1[i,k]<-log(p.abun[i,1]/p.abun[i,k])*p.abun[i,1]  ### first year is reference
   }}
   for (k in 2:nyears){
   D2[k]<-log(ttl_expected[k]/ttl_expected[1]) ### first year is reference
   D1[k]<--sum(d1[,k])
   D[k]<-D1[k]+D2[k]
   }
   prec <- 1/(sigma^2)
      sigma ~ dunif(0,6) 
  }  
  
    ",fill = TRUE)
sink()

# Bundle data
win.data.rf <- list(N = allbugsarray.mo.rf, nspp=dim(allbugsarray.mo.rf)[1], nmonths=dim(allbugsarray.mo.rf)[2], 
                 nyears=dim(allbugsarray.mo.rf)[3])

# Initial values
#inits<-function(){ list(logabun.pre = rep(1,18), logabun.post =rep(1,18))}
# Initial values
ntaxa.rf=53
nsite=1
nmonth=12
nyear=14
zst<-array(NA,dim=c(53,12,14))
for (m in 1:ntaxa.rf){
  for(i in 1:nmonth){
    for (k in 1:nyear){
      if(is.na(allbugsarray.mo[m,i,k]))
      {zst[m,i,k]<-NA
      next}
      zst[m,i,k]<-max(allbugsarray.mo[m,i,k],na.rm=T)
    }
  }
}

inits <- function(){ list(z = zst)}

# Parameters monitored

params <- c("D2", "D1", "D","expabun")

# MCMC settings
ni <- 12000
nt <- 4
nb <- 1200
nc <- 3

# Call JAGS from R (BRT 3 min)

out.bugs.ref1.rf<- jags(win.data.rf, inits, params, "coweetabugs_53_turnover.monthly.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
# Summarize posteriors
print(out.bugs.ref1.rf, dig = 2)

out.bugs.ref2.rf<- run.jags(model= "jags_turnover_model_monhtlybugs.R",
                         monitor=params,
                         data=win.data.rf,
                         inits=inits,
                         n.chains = nc, 
                         thin = nt, 
                         adapt=2e3,
                         sample = ni, 
                         burnin = nb,
                         summarise=F,
                         keep.jags.files = FALSE,
                         method = 'parallel'
)
out.bugs.summ.ref2.rf<-summary(out.bugs.ref2.rf)




bugstable.rf <- matrix(rep(NA,42*4), nrow=42)
bugstable.rf[1,2:4]<-0
bugstable.rf [15,2:4]<-0
bugstable.rf [29,2:4]<-0
bugstable.rf [2:14,2:4]<-out.bugs.summ.ref2.rf [1:13,c(4,1,3)]
bugstable.rf [16:28,2:4]<-out.bugs.summ.ref2.rf [14:26,c(4,1,3)]
bugstable.rf [30:42,2:4]<-out.bugs.summ.ref2.rf [27:39,c(4,1,3)]
bugstable.rf <-data.frame(bugstable.rf )
bugstable.rf [1:14,1]<-"D2"
bugstable.rf [15:28,1]<-"D1"
bugstable.rf [29:42,1]<-"D"
names(bugstable.rf )[1]<-"Dvar"
names(bugstable.rf )[2]<-"Mean"
names(bugstable.rf )[3]<-"L95"
names(bugstable.rf )[4]<-"U95"

year<-c(seq(1993,2006,1))
bugstable.rf<-cbind(bugstable.rf, year)
exp.abun.rf<-out.bugs.ref1.rf$BUGSoutput$mean$expabun
exp.abun.rf<-as.data.frame(exp.abun.rf)
exp.abun.rf.t<-t(exp.abun.rf)  ## expected abundances, years x spp

exp.abundance.rf<-exp.abun.rf.t
colnames(exp.abundance.rf)<-bugnames.rf[1,]
#write.csv(exp.abundance, file="exp_abund_53bc_annualmeans.csv")

exp.abundance.rf<-data.frame(exp.abundance.rf)
exp.abundance.rf<-cbind(exp.abundance.rf, year)
exp.abundance.rf$habitat<-"RF"
exp.abundance.rf<-exp.abundance.rf[,c(54,55,1:53)]


exp.abundance.all<-plyr::rbind.fill(exp.abundance,exp.abundance.rf)
#write.csv(exp.abundance.all, file="exp_abund_53bcrf_annualmeans_abundance.csv")





ggplot(bugstable.rf, aes(group=Dvar, color=Dvar))+
  #geom_point()+
  geom_line(data=bugstable.rf, aes(y=L95, x=year, group=Dvar), linetype=2)+
  geom_line(data=bugstable.rf, aes(y=U95, x=year, group=Dvar), linetype=2)+
  geom_line(data=bugstable.rf, aes(x=year, y=Mean, group=Dvar, color=Dvar,),size=1.25) +
  geom_point(data=bugstable.rf, aes(x=year, y=Mean, group=Dvar, color=Dvar,),size=3) +
  scale_color_manual(values=c("darkgreen","goldenrod2","dodgerblue3"))+
  scale_x_continuous(breaks=seq(1993,2006,1), limits=c(1993,2006))+
  theme_classic(base_size = 16)+
  theme(strip.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none",
        strip.text.x = element_text(size = 20),
        legend.background = element_blank(), 
        legend.title = element_blank())+
  geom_hline(yintercept = 0, linetype="dashed", size=1)+
  labs(y="Turnover", x="Year")

bugstable$habitat<-"BC"
bugstable.rf$habitat<-"RF"
bugstable.all<-rbind(bugstable, bugstable.rf)
#write.csv(bugstable.all, file="jagsturnovervalues_abundance_53bcrf.csv")

biomturn<-read.csv("jagsturnovervalues_53bcrf.csv")[-1]
turn.bay.rfbc<-ggplot(biomturn, aes(group=Dvar, color=Dvar))+
  #geom_point()+
  geom_line(data=biomturn, aes(y=L95, x=year, group=Dvar), linetype=2)+
  geom_line(data=biomturn, aes(y=U95, x=year, group=Dvar), linetype=2)+
  geom_line(data=biomturn, aes(x=year, y=Mean, group=Dvar, color=Dvar,),size=1.5) +
  geom_point(data=biomturn, aes(x=year, y=Mean, group=Dvar, color=Dvar,),size=3) +
  facet_wrap(.~habitat ,ncol=1,labeller = labeller(habitat = c("RF" = "Rockface ","BC" = "Mixed substrate")))+
  scale_color_manual(values=c("darkgreen","goldenrod2","dodgerblue3"),labels=c("D" = expression(italic("D")), "D1"=expression(italic(D[1])), "D2"=expression(italic(D[2]))))+
  scale_x_continuous(breaks=seq(1993,2006,1), limits=c(1993,2006))+
  scale_y_continuous(breaks=seq(-2, 1, .5), limits=c(-2,1))+
  theme_classic(base_size = 16)+
  theme(strip.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = c(.1,.5),
        strip.text.x = element_text(size = 20),
        legend.background = element_blank(), 
        legend.title = element_blank())+
  geom_hline(yintercept = 0, linetype="dashed", size=1)+
  labs(y="Turnover", x="Year")
ggsave("ws53_annualturnover_bcrf_bayesian_abundance_biggeryrange_25jul2025.jpeg", plot=turn.bay.rfbc, width=5, height=6)