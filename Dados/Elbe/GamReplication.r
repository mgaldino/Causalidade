# 
# replication code for distance plots in:
# How Lasting is Voter Gratitude? An Analysis of the Short- and Long-term Electoral Returns to Beneficial Policy"
# M. Bechtel and J. Hainmueller */

rm(list=ls())
library(foreign)
library(gam)


d <- read.dta("1998_2002.dta")
# reduce dataset to relevant variables
varnames <- c("wkr","bula","wkrname","Flooded","del_spd_z_vs","near_dist_allrivers",
              "ddpopdensity_","ddshforeign_","ddpopnetinp1000_","ddshpop_o60_","ddue_",
              "ddshagric_","ddshmanu_","ddshtradservice_","ddshotherservice_","ddsinc_SPD")
d <- na.omit(d[,varnames])
            
# fit GAM
gam.out     <- gam(del_spd_z_vs  ~s(near_dist_allrivers,4)+ddpopdensity_+ddshforeign_+ddpopnetinp1000_+ddshpop_o60_+ddue_
              +ddshagric_+ddshmanu_+ddshtradservice_+ddshotherservice_+ddsinc_SPD,data=d)

# left plot
plot(gam.out,terms="s(near_dist_allrivers, 4)",se=T,residuals=F,ylim=c(-12.5,11),
xlab="Distance to Elbe or Flooded Tributary (in km)",ylab="GAM Residual: Change in SPD PR Vote Share 2002-1998")
of   <- predict(gam.out,type="terms") 
pres <- residuals(gam.out)+of[,1]
points(y=pres,x=d$near_dist_allrivers,pch=20,col=gray(.7) )
points(y=pres[d$Flooded==1],x=d$near_dist_allrivers[d$Flooded==1],pch=20,col="black" )
legend("topright",legend=c("GAM Main-Effect Function","Pointwise CIs (2xSE)",
"Flooded","Other Districts"
),pch=c(NA,NA,20,20),lty=c(1,3,NA,NA),col=c("black","black","black",gray(.7)),cex=.75)

# right plot
plot(gam.out,terms="s(near_dist_allrivers, 4)",se=F,residuals=F,ylim=c(-12.5,11),
xlab="Distance to Elbe or Flooded Tributary (in km)",ylab="GAM Residual: Change in SPD PR Vote Share 2002-1998")
points(y=pres,x=d$near_dist_allrivers,pch=20,col=gray(.7))
pp <- c(0:3,5,6,8)
sname <- c("Schleswig-Holstein","Niedersachsen","Brandenburg","Mecklenburg-Vorpommern","Sachsen","Sachsen-Anhalt","ThŸringen")
for(i in 1:length(sname)){
yyy <- pres[d$bula==sname[i]]
xxx <- d$near_dist_allrivers[d$bula==sname[i]]
mod1 <- lm(yyy~xxx)
rrange <- seq(min(xxx),max(xxx),by=12)
locfit <- predict(mod1,newdata=data.frame(xxx=rrange))
points(rrange,locfit,col="black",pch=pp[i],type="o",cex=.75)
}
sname[7] <- "Thueringen"
leg <- paste(sname, c(rep("(West)",2),rep("(East)",5)))
legend("topright",legend=leg,col="black",pch=pp,lty=1,lwd=1,cex=.75)