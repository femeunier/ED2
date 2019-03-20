#==========================================================================================#
#==========================================================================================#
#     Leave these commands at the beginning.  They will refresh the session.               #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
#==========================================================================================#
#==========================================================================================#

here           = getwd()    # Current directory.
srcdir         = "/home/femeunier/Documents/ED2/R-utils" # Source  directory

there          = c("/home/femeunier/Documents/ED2/ED/run/analy",
                   "/home/femeunier/Documents/ED2/ED/run/analy/no_hydro",
                   "/home/femeunier/Documents/ED2/ED/run/analy/hydro1",
                   "/home/femeunier/Documents/ED2/ED/run/analy/hydro3")
                   # Directory where all analyses/history are 

monthbeg       = 1   # First month to useyearbeg        = thisyeara    # First year to consider
yearbeg        = 1500    # First year to consider
yearend        = 1500
places         = c("paracou","paracou","paracou","paracou")
names_simu     = c("No_liana","No_hydro","Hydro_1","Hydro_3")
LTY            = c(1,2,2,2)
colors         = c("black","green",'blue','red')
ntimes         = 6
slz.min        = -8.0 
sasmonth       =  seq(1,ntimes)

source(file.path(srcdir,"load.everything.r"))

datum_all=list()
for (i in seq(places)){

  print(paste('Reading',there[i]))
  place = places[i]
  thispoi = locations(where=place,here=there,yearbeg=yearbeg,yearend=yearend
                      ,monthbeg=monthbeg)
  inpref = file.path(there[i],place)
  
  datum      = create.monthly( ntimes  = ntimes
                               , montha  = monthbeg
                               , yeara   = yearbeg
                               , inpref  = inpref
                               , slz.min = slz.min) 
  datum = read.q.files(datum=datum,ntimes=ntimes,tresume=1,sasmonth=sasmonth)
  datum_all[[names_simu[i]]]=datum
}
plot.new()
par(mar=c(2,2,2,2),mfrow=c(1,2))

plot(x=0,y=0,xlim=c(1,ntimes),ylim=c(0,5))
ipft=18
for (i in seq(1,length(places))){
  lines(datum_all[[i]]$szpft$transp[,ndbh+1,ipft],col=colors[i],lty=LTY[i])
}

plot(x=0,y=0,xlim=c(1,ntimes),ylim=c(-2,2))

for (i in seq(2,length(places))){
  #ref=datum_all[[names_simu[1]]]$emean$gpp
  ref=datum_all[[1]]$szpft$transp[,ndbh+1,ipft]
  lines(ref-datum_all[[i]]$szpft$transp[,ndbh+1,ipft],col=colors[i],lty=LTY[i])
}

##########################################################################

plot.new()
par(mar=c(2,2,2,2),mfrow=c(1,2))

plot(x=0,y=0,xlim=c(1,ntimes),ylim=c(0,5))
ipft=17
for (i in seq(1,length(places))){
  lines(datum_all[[i]]$szpft$lai[,ndbh+1,npft+1],col=colors[i],lty=LTY[i])
}

plot(x=0,y=0,xlim=c(1,ntimes),ylim=c(0,2))

for (i in seq(1,length(places))){
  #ref=datum_all[[names_simu[1]]]$emean$gpp
  lines(datum_all[[i]]$szpft$lai[,ndbh+1,ipft],col=colors[i],lty=LTY[i])
}

##########################################################################

plot.new()
par(mar=c(2,2,2,2),mfrow=c(1,2))

plot(x=0,y=0,xlim=c(1,ntimes),ylim=c(0,15))
ipft=17
for (i in seq(1,length(places))){
  lines(datum_all[[i]]$szpft$agb[,ndbh+1,npft+1],col=colors[i],lty=LTY[i])
}

plot(x=0,y=0,xlim=c(1,ntimes),ylim=c(0,1))

for (i in seq(1,length(places))){
  #ref=datum_all[[names_simu[1]]]$emean$gpp
  lines(datum_all[[i]]$szpft$agb[,ndbh+1,ipft],col=colors[i],lty=LTY[i])
}
