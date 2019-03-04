#==========================================================================================#
#==========================================================================================#
#     Leave these commands at the beginning.  They will refresh the session.               #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
#==========================================================================================#
#==========================================================================================#

here           = getwd()    # Current directory.
there          = "/home/femeunier/Documents/ED2/ED/run"  # Directory where analyses/history are 
srcdir         = "/home/femeunier/Documents/ED2/R-utils" # Source  directory

monthbeg       = 1   # First month to useyearbeg        = thisyeara    # First year to consider
yearbeg        = 1500    # First year to consider
yearend        = 1500
place          = "paracou"
ntimes         = 15
slz.min        = -8.0 

sasmonth       =  seq(1,ntimes)

source(file.path(srcdir,"load.everything.r"))

thispoi = locations(where=place,here=there,yearbeg=yearbeg,yearend=yearend
                    ,monthbeg=monthbeg)
inpref = file.path(there,'analy',place)

datum      = create.monthly( ntimes  = ntimes
                             , montha  = monthbeg
                             , yeara   = yearbeg
                             , inpref  = inpref
                             , slz.min = slz.min
                             ) 
datum = read.q.files(datum=datum,ntimes=ntimes,tresume=1,sasmonth=sasmonth)

par(mar=c(6,4,4,4))
plot.new()

ipft=2
itime=ntimes
pft=datum$cohort$pft[[itime]]
N=length(pft)
N
plot(datum$cohort$dbh[[1]][datum$cohort$pft[[1]]==ipft],datum$cohort$height[[1]][datum$cohort$pft[[1]]==ipft],col='red')
lines(datum$cohort$dbh[[itime]][pft==ipft],datum$cohort$height[[itime]][pft==ipft],type='p')

