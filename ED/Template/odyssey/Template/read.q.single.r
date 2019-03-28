#==========================================================================================#
#==========================================================================================#
#     Leave these commands at the beginning.  They will refresh the session.               #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
#==========================================================================================#
#==========================================================================================#

here           = getwd()    # Current directory.
there          = "/home/femeunier/Documents/ED2/ED/run/analy"  # Directory where analyses/history are 
srcdir         = "/home/femeunier/Documents/ED2/R-utils" # Source  directory

monthbeg       = 1   # First month to useyearbeg        = thisyeara    # First year to consider
yearbeg        = 1500    # First year to consider
place          = "paracou"
ntimes         = 16
slz.min        = -8.0 

sasmonth       =  seq(1,ntimes)

source(file.path(srcdir,"load.everything.r"))

thispoi = locations(where=place,here=there,yearbeg=yearbeg,monthbeg=monthbeg)
inpref = file.path(there,place)

datum      = create.monthly( ntimes  = ntimes
                             , montha  = monthbeg
                             , yeara   = yearbeg
                             , inpref  = inpref
                             , slz.min = slz.min
                             ) 

datum = read.q.files(datum=datum,ntimes=ntimes,tresume=1,sasmonth=sasmonth)

ipft=17

# h5file <- '/home/femeunier/Documents/ED2/ED/run/histo/bci-S-2004-01-01-000000-g01.h5'
# 
# mymont    <- lapply(h5read_opt(h5file),FUN=aperm)
# names(mymont) <- gsub(x = names(mymont), pattern = "\\_", replacement = ".") 
# 
# ipaconow <- rep(sequence(mymont$NPATCHES.GLOBAL),times=mymont$PACO.N)
# PFT0 <- mymont$PFT
# HITE0 <- mymont$HITE
# DDBH0 <- mymont$DELTA.DBH
# DBH0 <- mymont$DBH
# plot(DBH0[PFT0==ipft],HITE0[PFT0==ipft],col='black',type='p')

itime=1
PFT=datum$cohort$pft[[itime]]
HITE=datum$cohort$height[[itime]]
DBH=datum$cohort$dbh[[itime]]
plot(DBH[PFT==ipft],HITE[PFT==ipft],ylim=c(0,36),col='black',type='p',xlab="DBH",ylab="HEIGHT",xlim=c(0,30))

itime=13
PFT2=datum$cohort$pft[[itime]]
HITE2=datum$cohort$height[[itime]]
DBH2=datum$cohort$dbh[[itime]]
lines(DBH2[PFT2==ipft],HITE2[PFT2==ipft],col='red',type='p')

itime=ntimes
PFT2=datum$cohort$pft[[itime]]
HITE2=datum$cohort$height[[itime]]
DBH2=datum$cohort$dbh[[itime]]
lines(DBH2[PFT2==ipft],HITE2[PFT2==ipft],col='green',type='p')

Npatches <- length(unique(datum$cohort$ipa[[1]]))
x0 <- c(0,0)
x1 <- c(0,0)

for (ipa in seq(Npatches)){
  pos1 <- which(datum$cohort$ipa[[1]] == ipa) 
  DBH1 <- datum$cohort$dbh[[1]][pos1][datum$cohort$pft[[1]][pos1]==17]
  HITE1 <- datum$cohort$height[[1]][pos1][datum$cohort$pft[[1]][pos1]==17]
  
  pos2 <- which(datum$cohort$ipa[[ntimes]] == ipa) 
  DBH2_temp <- datum$cohort$dbh[[ntimes]][pos2][datum$cohort$pft[[ntimes]][pos2]==17]
  HITE2_temp <- datum$cohort$height[[ntimes]][pos2][datum$cohort$pft[[ntimes]][pos2]==17]
  
  Ncohorts <- length(DBH1)
  
  if (Ncohorts>0){
    DBH2 <- DBH2_temp[1:(Ncohorts)]
    HITE2 <- HITE2_temp[1:(Ncohorts)]
    
    for (icohort in seq(Ncohorts)){
      dist <- colSums((rbind(rep(DBH1[icohort],Ncohorts),rep(HITE1[icohort],Ncohorts)) - rbind(DBH2,HITE2))^2)
      dist_min <- which(dist == min(dist))
      x0 <- rbind(x0,c(DBH1[icohort],HITE1[icohort]))
      x1 <- rbind(x1,c(DBH2[dist_min],HITE2[dist_min]))
    }
  }
}
arrows(x0 = x0[2:dim(x0)[1],1], y0 = x0[2:dim(x0)[1],2], x1 = x1[2:dim(x1)[1],1], y1 = x1[2:dim(x1)[1],2],
       length=0.05)
dbhs=seq(0,30,length.out = 1000)
b1ht=0.11
b2ht=0.87
hmax=35
href=61.7
h <- pmin(hmax,href*(1-exp(-b1ht*(dbhs**b2ht))))
lines(dbhs,h,lty=2,col='black')

b1ht=0.035
b2ht=0.69
h_tree <- pmin(hmax,href*(1-exp(-b1ht*(dbhs**b2ht))))
lines(dbhs,h_tree,lty=2,col='darkgreen')

##########################################################################

cohort <- datum$cohort

plot.new()
itime <- ntimes
ipft <- 17
pft <- cohort$pft[[itime]]
plot(cohort$dbh[[itime]][pft==ipft]-cohort$delta_DBH[[itime]][pft==ipft],cohort$height[[itime]][pft==ipft])
pft <- cohort$pft[[1]]
lines(cohort$dbh[[1]][pft==ipft]-cohort$delta_DBH[[1]][pft==ipft],cohort$height[[1]][pft==ipft],col='red',type='p')

dbhs=seq(0,30,length.out = 1000)
b1ht=0.11
b2ht=0.87
hmax=35
href=61.7
h <- pmin(hmax,href*(1-exp(-b1ht*(dbhs**b2ht))))
lines(dbhs,h,lty=2,col='black')


###########################################################################

itime=1
ipft=17
ipa=10

patch<-datum$patch
PA=datum$cohort$ipa[[itime]]
pos=(PA==ipa)
PFT=datum$cohort$pft[[itime]][pos]
HITE=datum$cohort$height[[itime]][pos]
DHITE=datum$cohort$delta_hite[[itime]][pos]
DBH=datum$cohort$dbh[[itime]][pos]

plot(DBH[PFT==ipft],HITE[PFT==ipft],ylim=c(0,36),xlim=c(0,200))

itime=ntimes
PA=datum$cohort$ipa[[itime]]
pos=(PA==ipa)
PFT2=datum$cohort$pft[[itime]][pos]
HITE2=datum$cohort$height[[itime]][pos]
DHITE=datum$cohort$delta_hite[[itime]][pos]
DBH2=datum$cohort$dbh[[itime]][pos]

lines(DBH2[PFT2==ipft],HITE2[PFT2==ipft],col='red',type='p')

################################################################

tapply(X=nplant[pos],INDEX = PFT[pos],FUN=sum)
sum(tapply(X=nplant[pos],INDEX = PFT[pos],FUN=sum))
tapply(X=LAI[pos],INDEX = PFT[pos],FUN=sum)
tapply(X=HITE[pos],INDEX = PFT[pos],FUN=mean)
tapply(X=PFT[pos],INDEX = PFT[pos],FUN=length)


###########################################################################

par(mar=c(6,4,4,4))
plot.new()

plot(datum$szpft$lai[,ndbh+1,2],type='l',ylim=c(0,6),col='lightgreen')
lines(datum$szpft$lai[,ndbh+1,3],type='l',ylim=c(0,6),col='green')
lines(datum$szpft$lai[,ndbh+1,4],type='l',ylim=c(0,6),col='darkgreen')
lines(datum$szpft$lai[,ndbh+1,17],type='l',ylim=c(0,6),col='blue')
lines(datum$szpft$lai[,ndbh+1,18],type='l',ylim=c(0,6),col='black')

plot(datum$szpft$agb[,ndbh+1,2],type='l',ylim=c(0,15),col='lightgreen')
lines(datum$szpft$agb[,ndbh+1,3],type='l',ylim=c(0,15),col='green')
lines(datum$szpft$agb[,ndbh+1,4],type='l',ylim=c(0,15),col='darkgreen')
lines(datum$szpft$agb[,ndbh+1,17],type='l',ylim=c(0,15),col='blue')
lines(datum$szpft$agb[,ndbh+1,18],type='l',ylim=c(0,15),col='black')

plot(rowSums(datum$szpft$nplant[,1:ndbh,2]),type='l',ylim=c(0,1),col='lightgreen')
lines(rowSums(datum$szpft$nplant[,1:ndbh,3]),type='l',col='green')
lines(rowSums(datum$szpft$nplant[,1:ndbh,4]),type='l',col='darkgreen')
lines(rowSums(datum$szpft$nplant[,1:ndbh,17]),type='l',col='blue')
lines(rowSums(datum$szpft$nplant[,1:ndbh,18]),type='l',col='black')

plot(datum$szpft$gpp[,ndbh+1,2],type='l',ylim=c(0,3),col='lightgreen')
lines(datum$szpft$gpp[,ndbh+1,3],type='l',col='green')
lines(datum$szpft$gpp[,ndbh+1,4],type='l',col='darkgreen')
lines(datum$szpft$gpp[,ndbh+1,17],type='l',col='blue')
lines(datum$szpft$gpp[,ndbh+1,18],type='l',col='black')

plot(datum$szpft$npp[,ndbh+1,2],type='l',ylim=c(0,1),col='lightgreen')
lines(datum$szpft$npp[,ndbh+1,3],type='l',col='green')
lines(datum$szpft$npp[,ndbh+1,4],type='l',col='darkgreen')
lines(datum$szpft$npp[,ndbh+1,17],type='l',col='blue')
lines(datum$szpft$npp[,ndbh+1,18],type='l',col='black')

plot(datum$szpft$zRWU[,ndbh+1,2],type='l',ylim=c(slz.min,0),col='lightgreen')
lines(datum$szpft$zRWU[,ndbh+1,3],type='l',col='green')
lines(datum$szpft$zRWU[,ndbh+1,4],type='l',col='darkgreen')
lines(datum$szpft$zRWU[,ndbh+1,17],type='l',col='blue')
lines(datum$szpft$zRWU[,ndbh+1,18],type='l',col='black')


plot(datum$emean$agb/datum$emean$agb[1])

plot(datum$szpft$lai[,ndbh+1,4])

plot(datum$cohort$dbh[[1]],datum$cohort$bstorage[[1]])
lines(datum$cohort$dbh[[6]],datum$cohort$bstorage[[6]],type='p',col='red')

datum$szpft$lai[,ndbh+1,2]

plot.new()
ipft=17
pft=datum$cohort$pft[[itime]]
N=length(pft)
N
plot(datum$cohort$dbh[[1]][datum$cohort$pft[[1]]==ipft],datum$cohort$height[[1]][datum$cohort$pft[[1]]==ipft],col='red')
lines(datum$cohort$dbh[[itime]][pft==ipft],datum$cohort$height[[itime]][pft==ipft],type='p')

itime=ntimes
Npatches=max(datum$cohort$ipa[[itime]])
maxH_L=maxH_T=matrix(NA,Npatches)
for (i in seq(1,Npatches)){
  h=datum$cohort$height[[itime]][datum$cohort$ipa[[itime]]==i]
  pft=datum$cohort$pft[[itime]][datum$cohort$ipa[[itime]]==i]
  if (any(pft==17)){
    maxH_L[i]=max(h[pft==17])
  }
  if (any(pft!=17)){
  maxH_T[i]=max(h[pft!=17])
  }
}


plot(seq(1,Npatches),maxH_T,type='p',col='green')
lines(seq(1,Npatches),maxH_L,type='p',col='blue')

b1Rd  = -1.1140580
b2Rd  =  0.4223014

ipatch=36
datum$patch$transp[[itime]][ipatch]
pos=datum$cohort$ipa[[itime]]==ipatch
h=datum$cohort$height[[itime]][pos]
pft=datum$cohort$pft[[itime]][pos]
lai=datum$cohort$lai[[itime]][pos]
dbh=datum$cohort$dbh[[itime]][pos]
n=datum$cohort$nplant[[itime]][pos]
broot=datum$cohort$broot[[itime]][pos]

hite=h[pft==17]
root_depth = b1Rd * hite ** b2Rd

SLZ     =   c(-8.000,-6.959,-5.995,-5.108,-4.296,-3.560,-2.897,-2.307,-1.789,-1.340,-0.961,-0.648,-0.400,-0.200,-0.100,-0.050)
nzg =length(SLZ)
slz8=SLZ
root_beta   =   0.001
krdepth=12
root_frac=matrix(0.,16)


for (k in seq(krdepth,16)){
  current_layer_depth = -slz8[k]
  if (k+1 <= nzg){
    above_layer_depth = -slz8[k+1]
  }  else {
    above_layer_depth = 0.
  }
  root_frac[k] = ((root_beta ** (above_layer_depth   / (-slz8[krdepth]))
                -(root_beta) ** (current_layer_depth / (-slz8[krdepth]))))
}
plot(dbh[pft==17],h[pft==17])