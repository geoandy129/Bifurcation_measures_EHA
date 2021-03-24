## Bifurcation measurements in EHA
##
## Andy Mitten, Keele University
## a.j.mitten@keele.ac.uk
#
# This code enables you to measure the bifurcations in an evolutive harmonic
# analysis (EHA) by selecting a target location within the dataset being 
# analysed and determine the normalised change in amplitude from the MTM
# periodograms used to generate the EHA, thus quantifying hiatus duration.
#
# The code consists of the following steps:
#
# 1. Load in the data by setting up the working directory and data set. The
#    the dataset should be in a .txt tab deliminated format and should not 
#    contain headings. The first column of the data should be locational
#    (hieght/depth/time), and the second should be the dataset.
# 2. Basic plotting - Generate the MTM periodogram and confidence intervals for
#    the entire section. This section also enables you to generate an evolutive 
#    harmoinic analysis for the dataset. I suggest saving this as a .pdf file
#    as it may be useful during the picking process.
# 3. Identify the target interval (where the bifurcation is located within the 
#    dataset). Set the top and basal locations.
# 4. Identify the maximum and minimum frequencies you would like to search for
#    birfuraction peaks within. The minimum should be lower and the maximum 
#    should be higher than the bifurcation, respectively.
# 5. Check you have the right sampled frequencies and the correct target depth
#    by running a second "high resolution" EHA. Use this EHA to identify target
#    and error frequencies.
# 6. Run the code line by line carefully where instructed and the carefully read
#    the instructions below for the picking of the amplitudes.
# 7. Complete the amplitude difference matrix and calculate your hiatus 
#    duration and error duration. This will also give you the error associated 
#    in that pick as a percentage.
#
# The base code used for this and the package it is built within is Astrochron 
# for R (Myers, 2014). This code is simply using aspects of Astrochron and 
# streamlining the process of bifurcation picking. Full credit a citation should
# be as follows:
# Meyers, S.R. (2014). Astrochron: An R Package for Astrochronology. 
# https://cran.r-project.org/package=astrochron

rm(list = ls())


#loading------------------------------------------------------------------------

library(astrochron)
#SET WORKING DIRECTORY
D=setwd("WHATEVER THIS IS FOR YOU")
#CHANGE FILE NAME TO WHAT YOU WNAT TO OBSERVE
df.Faulted=read(file="YOUR DATA FILE.txt",d=0,h="no",skip=0,srt=T,ave=T,check=T,
                genplot=T,verbose=T)

#Basic plotting-----------------------------------------------------------------

pdf(file = "THE LOCATION AND NAME YOU WANT TO SAVE AS.pdf")
a=mtm(df.Faulted,tbw=2,ar=FALSE,pl=2,detrend=T,demean=T,
      xmin=0,xmax=0.5,genplot=T,sigID=F,verbose=F)
dev.off()

b=eha(df.Faulted,tbw=2,win=30,pad=10000,output=3,fmin=0,fmax=0.5,
      detrend=T,demean=T,step=1,palette=1,ncolors=2000,
      verbose=F,xlab="Frequency",ylab="Time (kyr)",
      genplot=2,siglevel=0.9,centerZero=T,sigID=F)


# Setting up the search---------------------------------------------------------

minL=40 #Basal location of bifurcation
maxL=60 #Top location of bifurcation
itL=(maxL-minL)/20 #DO NOT CHANGE THIS

# cHECK YOUR TARGET FREQUENCY---------------------------------------------------

toSearch=seq(minL,maxL,itL) #DO NOT CHANGE THIS
intrestFreqmin=0.12 #minimum search frequency
intrestFreqmax=0.28 #maximum search frequency

yMAx=2 #Maximum values from y-axis in MTM + a padding... see what works!

c=eha(df.Faulted,tbw=2,win=30,pad=10000,output=3,fmin=intrestFreqmin
      ,fmax=intrestFreqmax,detrend=T,demean=T,step=1,
      palette=1,ncolors=2000,
      verbose=F,xlab="Frequency",ylab="Time (kyr)",
      genplot=2,siglevel=0.9,centerZero=T,sigID=F)

# Adjust the target frequency to read the most statistically significant of the 
# bifurcation

signifFreq=0.14 # Used to estimate possible error
originalFreq=0.2 # Used to determine hiatus duration

# DOES EVERYTHING LOOK RIGHT? IF NOT CHNAGE THE VALUES ABAOVE TO MAKE SURE 
# YOU HIT YOUR TARGET.
#
# SAVE THIS AS A PDF.
#
# RUN THIS CAREFULLY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# READ THIS IN DETAIL.
# 
# BEFORE GOING ANY FURTHER
# DO THE FOLLOWING:
# 1- CLOSE ALL EXTERNAL WINDOWS, KEEP THOSE IN THE PLOT PANE
# 2- ENSURE YOU CAN SEE THE SECOND "ZOOMED-IN" EHA IN A PDF READER.
# 3- CLICK CNTRL+ENTR 3 TIMES THIS WILL TAKE YOU FROM SPECn TO AMPn
# 4- SELECT YOUR POINTS
# 5- CLICK FINISH IN THE TOP RIGHT OF THE PLOTS PANE
# 6- REPEAT #3-#5 UNTIL YOU REACH THE END OF THE SECTION.

# RUN THIS CAREFULLY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# MEASURING---------------------------------------------------------------------

spec1=extract(b,get=toSearch[1],genplot=F)
freqs1=peak(spec1,genplot=F)
amp1=idPts(freqs1,ptsize=5,xmin=intrestFreqmin,
           xmax=intrestFreqmax,ymin=0,ymax=yMAx,
      logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

spec2=extract(b,get=toSearch[2],genplot=F)
freqs2=peak(spec2,genplot=F)
amp2=idPts(freqs2,ptsize=5,xmin=intrestFreqmin,
           xmax=intrestFreqmax,ymin=0,ymax=yMAx,
           logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

spec3=extract(b,get=toSearch[3],genplot=F)
freqs3=peak(spec3,genplot=F)
amp3=idPts(freqs3,ptsize=5,xmin=intrestFreqmin,
           xmax=intrestFreqmax,ymin=0,ymax=yMAx,
           logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

spec4=extract(b,toSearch[4],genplot=F)
freqs4=peak(spec4,genplot=F)
amp4=idPts(freqs4,ptsize=5,xmin=intrestFreqmin,
           xmax=intrestFreqmax,ymin=0,ymax=yMAx,
           logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

spec5=extract(b,get=toSearch[5],genplot=F)
freqs5=peak(spec5,genplot=F)
amp5=idPts(freqs5,ptsize=5,xmin=intrestFreqmin,
           xmax=intrestFreqmax,ymin=0,ymax=yMAx,
           logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

spec6=extract(b,get=toSearch[6],genplot=F)
freqs6=peak(spec6,genplot=F)
amp6=idPts(freqs6,ptsize=5,xmin=intrestFreqmin,
           xmax=intrestFreqmax,ymin=0,ymax=yMAx,
           logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

spec7=extract(b,get=toSearch[7],genplot=F)
freqs7=peak(spec7,genplot=F)
amp7=idPts(freqs7,ptsize=5,xmin=intrestFreqmin,
           xmax=intrestFreqmax,ymin=0,ymax=yMAx,
           logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

spec8=extract(b,get=toSearch[8],genplot=F)
freqs8=peak(spec8,genplot=F)
amp8=idPts(freqs8,ptsize=5,xmin=intrestFreqmin,
           xmax=intrestFreqmax,ymin=0,ymax=yMAx,
           logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

spec9=extract(b,get=toSearch[9],genplot=F)
freqs9=peak(spec9,genplot=F)
amp9=idPts(freqs9,ptsize=5,xmin=intrestFreqmin,
           xmax=intrestFreqmax,ymin=0,ymax=yMAx,
           logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

spec10=extract(b,get=toSearch[10],genplot=F)
freqs10=peak(spec10,genplot=F)
amp10=idPts(freqs10,ptsize=5,xmin=intrestFreqmin,
           xmax=intrestFreqmax,ymin=0,ymax=yMAx,
           logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

spec11=extract(b,get=toSearch[11],genplot=F)
freqs11=peak(spec11,genplot=F)
amp11=idPts(freqs11,ptsize=5,xmin=intrestFreqmin,
            xmax=intrestFreqmax,ymin=0,ymax=yMAx,
            logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

spec12=extract(b,get=toSearch[12],genplot=F)
freqs12=peak(spec12,genplot=F)
amp12=idPts(freqs12,ptsize=5,xmin=intrestFreqmin,
           xmax=intrestFreqmax,ymin=0,ymax=yMAx,
           logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

spec13=extract(b,get=toSearch[13],genplot=F)
freqs13=peak(spec13,genplot=F)
amp13=idPts(freqs13,ptsize=5,xmin=intrestFreqmin,
           xmax=intrestFreqmax,ymin=0,ymax=yMAx,
           logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

spec14=extract(b,get=toSearch[14],genplot=F)
freqs14=peak(spec14,genplot=F)
amp14=idPts(freqs14,ptsize=5,xmin=intrestFreqmin,
           xmax=intrestFreqmax,ymin=0,ymax=yMAx,
           logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

spec15=extract(b,get=toSearch[15],genplot=F)
freqs15=peak(spec15,genplot=F)
amp15=idPts(freqs15,ptsize=5,xmin=intrestFreqmin,
           xmax=intrestFreqmax,ymin=0,ymax=yMAx,
           logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

spec16=extract(b,get=toSearch[16],genplot=F)
freqs16=peak(spec16,genplot=F)
amp16=idPts(freqs16,ptsize=5,xmin=intrestFreqmin,
           xmax=intrestFreqmax,ymin=0,ymax=yMAx,
           logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

spec17=extract(b,get=toSearch[17],genplot=F)
freqs17=peak(spec17,genplot=F)
amp17=idPts(freqs17,ptsize=5,xmin=intrestFreqmin,
           xmax=intrestFreqmax,ymin=0,ymax=yMAx,
           logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

spec18=extract(b,get=toSearch[18],genplot=F)
freqs18=peak(spec18,genplot=F)
amp18=idPts(freqs18,ptsize=5,xmin=intrestFreqmin,
           xmax=intrestFreqmax,ymin=0,ymax=yMAx,
           logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

spec19=extract(b,get=toSearch[19],genplot=F)
freqs19=peak(spec19,genplot=F)
amp19=idPts(freqs19,ptsize=5,xmin=intrestFreqmin,
           xmax=intrestFreqmax,ymin=0,ymax=yMAx,
           logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

spec20=extract(b,get=toSearch[20],genplot=F)
freqs20=peak(spec20,genplot=F)
amp20=idPts(freqs20,ptsize=5,xmin=intrestFreqmin,
           xmax=intrestFreqmax,ymin=0,ymax=yMAx,
           logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

spec21=extract(b,get=toSearch[21],genplot=F)
freqs21=peak(spec21,genplot=F)
amp21=idPts(freqs21,ptsize=5,xmin=intrestFreqmin,
            xmax=intrestFreqmax,ymin=0,ymax=yMAx,
            logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)


# Fault throw calculations------------------------------------------------------

#DELTA MU CALCULATIONS

amp1D=(amp1[1,2]-amp1[2,2])/(amp1[1,2]+amp1[2,2])
amp2D=(amp2[1,2]-amp2[2,2])/(amp2[1,2]+amp2[2,2])
amp3D=(amp3[1,2]-amp3[2,2])/(amp3[1,2]+amp3[2,2])
amp4D=(amp4[1,2]-amp4[2,2])/(amp4[1,2]+amp4[2,2])
amp5D=(amp5[1,2]-amp5[2,2])/(amp5[1,2]+amp5[2,2])
amp6D=(amp6[1,2]-amp6[2,2])/(amp6[1,2]+amp6[2,2])
amp7D=(amp7[1,2]-amp7[2,2])/(amp7[1,2]+amp7[2,2])
amp8D=(amp8[1,2]-amp8[2,2])/(amp8[1,2]+amp8[2,2])
amp9D=(amp9[1,2]-amp9[2,2])/(amp9[1,2]+amp9[2,2])
amp10D=(amp10[1,2]-amp10[2,2])/(amp10[1,2]+amp10[2,2])
amp11D=(amp11[1,2]-amp11[2,2])/(amp11[1,2]+amp11[2,2])
amp12D=(amp12[1,2]-amp12[2,2])/(amp12[1,2]+amp12[2,2])
amp13D=(amp13[1,2]-amp13[2,2])/(amp13[1,2]+amp13[2,2])
amp14D=(amp14[1,2]-amp14[2,2])/(amp14[1,2]+amp14[2,2])
amp15D=(amp15[1,2]-amp15[2,2])/(amp15[1,2]+amp15[2,2])
amp16D=(amp16[1,2]-amp16[2,2])/(amp16[1,2]+amp16[2,2])
amp17D=(amp17[1,2]-amp17[2,2])/(amp17[1,2]+amp17[2,2])
amp18D=(amp18[1,2]-amp18[2,2])/(amp18[1,2]+amp18[2,2])
amp19D=(amp19[1,2]-amp19[2,2])/(amp19[1,2]+amp19[2,2])
amp20D=(amp20[1,2]-amp20[2,2])/(amp20[1,2]+amp20[2,2])
amp21D=(amp21[1,2]-amp21[2,2])/(amp21[1,2]+amp21[2,2])

AD=c(amp1D, amp2D,amp3D,amp4D,amp5D,amp6D,amp7D,amp8D,
     amp9D,amp10D,amp11D,amp12D,amp13D,amp14D,amp15D,amp16D,
     amp17D,amp18D,amp19D,amp20D,amp21D)

plot(AD,toSearch)

AD=data.frame(AD,toSearch)

# PICKING DELTA MU FOR YOUR DATASET---------------------------------------------

ADMax=idPts(AD,ptsize=3,xmin=-1,
            xmax=2,ymin=minL,ymax=maxL,
            logx=F,logy=F,plotype=2,annotate=1,output=1,verbose=T)

FMax=(0.826*ADMax)+0.494
print(FMax)

HiatusDuration=(1/originalFreq)*FMax[1,1]

ErrorHiatus=(1/signifFreq)*FMax[1,1]

print(HiatusDuration)
print(ErrorHiatus)

Error=100-((ErrorHiatus/HiatusDuration)*100)
print(Error)



