
#### as an example, counts of 53 species at the downstream-most 4 Etowah sites ("Yellow Creek 7, 8, 9, 10), 19 annual surveys
yc.7.10.counts <- read.csv("yc.7.10.counts.csv")
nspp<-53 ## including everybody -4 with 0 obs at these sites; 4 sites, 19 years

### need a 3 dim matrix, sites, spp, years
str(yc.7.10.counts)
yc.7.10.counts$site<-as.factor(yc.7.10.counts$site)
library(reshape2)
library(tidyverse)
allcounts<-melt(yc.7.10.counts)
str(allcounts)

allcountsarray<-array(NA, dim=c(53,4,19))
str(allcountsarray)

allcountsarray[,,]<-as.array(allcounts$value)
str(allcountsarray) 

###ta-da### check a couple of values - say counts site 9, P. nigo, 2008 should be 34
allcountsarray[46,3,12] #34  check
allcountsarray[,1,12] #all counts, 2008 , site 7 - these align; all good


library(rjags)
library(R2jags)
library(runjags)
library(tidyverse)
bugs.ann<-read.csv("bugs_annual.csv")
bugs.ann.bc<-subset(bugs.ann, stream=="55" &habitat=="BC")
drop.taxa<-c("scrapers", "shredders", "tot.invert", "inv.pred", "gatherers", "filterers", "tot.chordata", "desmog", "eurycea","inv.predB", "cambarus")
bugs.ann.bc<-bugs.ann.bc[ , -which(names(bugs.ann.bc) %in% drop.taxa)]


bugs.ann.bc53<-subset(bugs.ann, stream=="53" &habitat=="BC")
drop.taxa<-c("scrapers", "shredders", "tot.invert", "inv.pred", "gatherers", "filterers", "tot.chordata", "desmog", "eurycea","inv.predB", "cambarus")
bugs.ann.bc53<-bugs.ann.bc53[ , -which(names(bugs.ann.bc53) %in% drop.taxa)]
bugs.ann.l53<-bugs.ann.bc53%>%gather(taxa,totbiom,c(40:106))
bugs.ann.l53<-rename(bugs.ann.l53, taxon=taxa)

######Calculate turnover using abundance########

#####Monthly Abundance Data#####
abund.mo<-read.csv("abund.w.spaggregated.csv")
abund.mo<-abund.mo[1:672,]#data is duplicated for some reason
abund.taxa.53bc=subset(abund.mo, stream=="53"& habitat=="BC")
abund.taxa.53a<-abund.taxa.53bc[,c(15,17,30,32,39:105)]
abund.taxa.53a<-abund.taxa.53a[ , -which(names(abund.taxa.53a) %in% drop.taxa)]



month.order<-data.frame(month=c(1,2,3,4,5,6,7,8,9,10,11,12), sample.month=c("e.1","f.2","g.3","h.4","i.5","j.6","k.7", "l.8","a.9","b.10","c.11","d.12"))
abund.taxa.53b<-left_join(abund.taxa.53a,month.order, by="month")
abund.taxa.53b<-abund.taxa.53b[,c(1,2,3,67,4,5:65)]




abund.taxa.53l<-abund.taxa.53b%>%gather(taxa,totabund,c(6:66))
abund.taxa.53l<-rename(abund.taxa.53l, taxon=taxa)
#bc.taxa<-data.frame(unique(abund.taxa.53l$taxa))
#abund.taxa.53l.bc<-filter(abund.taxa.53l, habitat=="BC")
#abund.taxa.53l.rf<-filter(abund.taxa.53l, habitat=="RF")
#write.csv(file="bctaxa.sppaggregated.csv", bc.taxa)

allbugs<-bugs.ann.l53[,c(3,40:41)]
#allbugs.mo<-abund.taxa.53l[,c(2,4:6)]
allbugs.mo<-abund.taxa.53l[,c(6,4,5,7)]
allbugs.mo.rf<-abund.taxa.53l.rf[,c(6,4,5,7)]
#allbugs.mo$month<-as.factor(allbugs.mo$month)

#allbugs$sample.yr<-as.factor(allbugs$sample.yr)
#as.factor(allbugs$site<-1)
allbugs<-allbugs[,c(2,1,3)]

allbugs2<-allbugs%>%spread(sample.yr, totabund)
allbugs.mo2<-allbugs.mo%>%spread(sample.yr, totabund)
allbugs.mo2.rf<-allbugs.mo.rf%>%spread(sample.yr, totabund)

allbugs3<-melt(allbugs2)
allbugs.mo3<-melt(allbugs.mo2)
allbugs.mo3.rf<-melt(allbugs.mo2.rf)
allbugs.mo4<-allbugs.mo3[order(  allbugs.mo3$variable,allbugs.mo3$sample.month, allbugs.mo3$taxon),]
allbugs.mo4.rf<-allbugs.mo3.rf[order(  allbugs.mo3.rf$variable,allbugs.mo3.rf$sample.month, allbugs.mo3.rf$taxon),]
write.csv(allbugs.mo4, "53bc_inputjagsdataformatted_abundance.csv")

nbugs<-length(unique(allbugs3$taxon))
nbugs<-length(unique(allbugs.mo$taxon))
nyears<-length(unique(allbugs.mo$sample.yr))
no.months<-length(unique(allbugs.mo$sample.month))


allbugsarray<-matrix(NA, 61,14 )
allbugsarray.mo<-array(NA, dim=c(61,12,14))
str(allbugsarray)
str(allbugsarray.mo)

allbugsarray[,]<-as.matrix(allbugs3$value)
allbugsarray.mo[,,]<-as.array(allbugs.mo4$value)
str(allbugsarray.mo) 
bugnames<-matrix(NA, 1,61 )
str(bugnames)
bugnames[,]<-as.matrix(unique(allbugs.mo4$taxon))

###ta-da### check a couple of values 
allbugsarray.mo[9,1,1] #should be 75; corduleg, a.9, 1993
allbugsarray.mo[40,5,1] # should be 12.5; parapsyche, e.1, 1993






#### estimate community turnover as in Shimadzu et al.
#### let the species-specifc counts each year at the four sites be draws from an underlying mean 
#### expected abundance per shoal sample; could add a random effect for site?

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




#exp.abun.bc<-out.bugs.summ[-c(1:39),4]
#exp.abun<-as.data.frame(exp.abun.bc)
#exp.abun.bc2 <- as.data.frame(matrix(exp.abun[,1], byrow=TRUE, ncol = 14))
#exp.abun.bc.t<-t(exp.abun.bc2)  ## expected abundances, years x spp

#exp.abundance<-exp.abun.bc.t
#colnames(exp.abundance)<-bugnames[1,]
#exp.abundance<-data.frame(exp.abundance)
#exp.abundance<-cbind(exp.abundance, year)
#exp.abundance$habitat<-"BC"
#exp.abundance<-exp.abundance[,c(62,63,1:61)]

write.csv(exp.abundance, file="exp_abund_53bc_abundance_annualmeans.csv")




# compare D values to Shimadzu code

DD<-function(x,ref.t=1, zero.rm=FALSE){
  lmb<-apply(x,1,sum)
  D2<- log(lmb/lmb[ref.t])
  
  x.p<-t(apply(x,1,function(z)z/sum(z)))
  Pt<-x.p[ref.t,]
  if(zero.rm==FALSE){
    D1<- -t(apply(x.p,1, function(z)ifelse(Pt==0, 0, log(Pt/z)))) %*% Pt
  }else{
    D1<- -t(apply(x.p,1, function(z)ifelse(Pt==0|z==0, 0, log(Pt/z)))) %*% Pt
  }
  D<-D1 + D2
  data.frame(D, D1, D2)
}

DD(exp.abun.t, ref.t=1, zero.rm = FALSE)  ### D estimates using Shimadzu code
D.data  ### D estimates from jags code
### these are the same for D2, really  for D1 - maybe check those calculations

library(ggplot2)
ggplot(D.data, aes(x=year, y=D))+
  geom_point()+
  geom_line(color="red", size=1.25) +
  theme_classic(base_size = 12)+
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")+
  geom_line(aes(x=year, y=D1), color="green",size=1.25)+
  geom_line(aes(x=year, y=D2), color="blue",size=1.25)+
  geom_hline(yintercept = 0, linetype="dashed")+
  ylab("Community turnover")


ggplot(bugstable, aes(group=Dvar, color=Dvar))+
  #geom_point()+
  geom_line(data=bugstable, aes(y=L95, x=year, group=Dvar), linetype=2)+
  geom_line(data=bugstable, aes(y=U95, x=year, group=Dvar), linetype=2)+
  geom_line(data=bugstable, aes(x=year, y=Mean, group=Dvar, color=Dvar,),size=1.25) +
  geom_point(data=bugstable, aes(x=year, y=Mean, group=Dvar, color=Dvar,),size=3) +
  scale_color_manual(values=c("darkgreen","goldenrod2","dodgerblue3"))+
  scale_x_continuous(breaks=seq(1993,2006,1), limits=c(1993,2006))+
  theme_classic(base_size = 16)+
  theme(strip.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = c(.1,.9),
        strip.text.x = element_text(size = 20),
        legend.background = element_blank(), 
        legend.title = element_blank())+
  geom_hline(yintercept = 0, linetype="dashed", size=1)+
  labs(y="Turnover", x="Year")

ggplot(D.data, aes(x=year))+
  geom_line(aes(y=D1), color="blue") +
  geom_line(aes(y=D2)) +
  #geom_line(aes(y=D), color="darkgrey", linetype="twodash")+
  theme_classic(base_size = 12)+
  theme(strip.background = element_blank()) +
  geom_hline(yintercept = 0, linetype="dashed")+
  ylab("Community turnover")




jm.sum=summary(out.bugs.ref1)
jm.mcmc=as.mcmc(out.bugs.ref1)

nvalues<-1
year2<-c(seq(1994,2006,1))
#jm1_mcmc_combi<-as.mcmc(rbind(jm1.mcmc[[1]],jm1.mcmc[[2]], jm1.mcmc[[3]], jm1.mcmc[[4]], jm1.mcmc[[5]]  ))


plot(year2, jm.sum[65:77,"Mean"])


D2<-D.data %>%
  mutate(delD = D - lag(D, default = 0), delD1=D1 - lag(D1, default = 0),delD2=D2 - lag(D2, default = 0) )











######Rockface######
abund.taxa.53rf=subset(abund.mo, stream=="53"& habitat=="RF")
abund.taxa.53rfa<-abund.taxa.53rf[,c(15,17,30,32,39:105)]
abund.taxa.53rfa<-abund.taxa.53rfa[ , -which(names(abund.taxa.53rfa) %in% drop.taxa)]
abund.taxa.53rfa.nosum<-abund.taxa.53rfa[,c(5:66)]
abund.taxa.53rfa.nosum = abund.taxa.53rfa.nosum[,colSums(abund.taxa.53rfa.nosum) == 0]
abund.taxa.53rfa<-abund.taxa.53rfa[ , -which(names(abund.taxa.53rfa) %in% colnames(abund.taxa.53rfa.nosum))]



month.order<-data.frame(month=c(1,2,3,4,5,6,7,8,9,10,11,12), sample.month=c("e.1","f.2","g.3","h.4","i.5","j.6","k.7", "l.8","a.9","b.10","c.11","d.12"))
abund.taxa.53rfb<-left_join(abund.taxa.53rfa,month.order, by="month")
abund.taxa.53rfb<-abund.taxa.53rfb[,c(1,2,3,58,4,5:57)]




abund.taxa.53l.rf<-abund.taxa.53rfb%>%gather(taxa,totabund,c(6:58))
abund.taxa.53l.rf<-rename(abund.taxa.53l.rf, taxon=taxa)
rf.taxa<-data.frame(unique(abund.taxa.53l.rf$taxa))
write.csv(file="rftaxa.sppaggregated.csv", rf.taxa)

#abund.taxa.53l.bc<-filter(abund.taxa.53l, habitat=="BC")
#abund.taxa.53l.rf<-filter(abund.taxa.53l, habitat=="RF")




allbugs.mo.rf<-abund.taxa.53l.rf[,c(6,4,5,7)]
#allbugs.mo$month<-as.factor(allbugs.mo$month)



allbugs.mo2.rf<-allbugs.mo.rf%>%spread(sample.yr, totabund)


allbugs.mo3.rf<-melt(allbugs.mo2.rf)

allbugs.mo4.rf<-allbugs.mo3.rf[order(  allbugs.mo3.rf$variable,allbugs.mo3.rf$sample.month, allbugs.mo3.rf$taxon),]
write.csv(allbugs.mo4.rf, "53rf_inputjagsdataformatted_abundance.csv")

nbugs.rf<-length(unique(allbugs.mo.rf$taxon))
nyears.rf<-length(unique(allbugs.mo.rf$sample.yr))
no.months.rf<-length(unique(allbugs.mo.rf$sample.month))


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



#exp.abun.rf<-out.bugs.summ.ref2.rf[-c(1:39),4]
#exp.abun.rf<-as.data.frame(exp.abun.rf)
#exp.abun.rf2 <- as.data.frame(matrix(exp.abun.rf[,1], byrow=TRUE, ncol = 14))
#exp.abun.rf.t<-t(exp.abun.rf2)  ## expected abundances, years x spp

#exp.abundance.rf<-exp.abun.rf.t
#colnames(exp.abundance.rf)<-bugnames.rf[1,]
#exp.abundance.rf<-data.frame(exp.abundance.rf)
#exp.abundance.rf<-cbind(exp.abundance.rf, year)
#exp.abundance.rf$habitat<-"RF"
#exp.abundance.rf<-exp.abundance.rf[,c(54,55,1:53)]

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









write.csv(exp.abundance.rf, file="exp_abund_53rf_annualmeans.csv")


exp.abundance.all<-plyr::rbind.fill(exp.abundance,exp.abundance.rf)
write.csv(exp.abundance.all, file="exp_abund_53bcrf_annualmeans_abundance.csv")





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
write.csv(bugstable.all, file="jagsturnovervalues_abundance_53bcrf.csv")
biomturn<-read.csv("jagsturnovervalues_53bcrf.csv")[-1]

bugstable.all.biom.abund<-cbind(biomturn,bugstable.all)
names(bugstable.all.biom.abund )[7]<-"Dvar.a"
names(bugstable.all.biom.abund )[8]<-"Mean.a"
names(bugstable.all.biom.abund )[9]<-"L95.a"
names(bugstable.all.biom.abund )[10]<-"U95.a"
bugstable.all.biom.abund<-bugstable.all.biom.abund[,-c(11:12)]




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


ggplot(bugstable.all.biom.abund, aes(x=Mean, y=Mean.a,color=Dvar, group=Dvar))+
  geom_point(size=3)+
  scale_color_manual(values=c("darkgreen","goldenrod2","dodgerblue3"))+
  facet_wrap(.~habitat ,ncol=1,labeller = labeller(habitat = c("RF" = "Rock face ","BC" = "Mixed cobble")))
  
  
