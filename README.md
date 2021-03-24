# Bifurcation_measures_EHA

## Bifurcation measurements in EHA
##
## Andy Mitten, Keele University
## a.j.mitten@keele.ac.uk
#
# This code enables you to measure the bifurcations in an evolutive harmonic analysis (EHA) by selecting a target location within the dataset being analysed and determine the normalised change in amplitude from the MTM periodograms used to generate the EHA, thus quantifying hiatus duration.
#
# The code consists of the following steps:
#
# 1. Load in the data by setting up the working directory and data set. The the dataset should be in a .txt tab deliminated format and should not contain headings. The first column of the data should be locational (hieght/depth/time), and the second should be the dataset.
# 2. Basic plotting - Generate the MTM periodogram and confidence intervals for the entire section. This section also enables you to generate an evolutive harmoinic analysis for the dataset. I suggest saving this as a .pdf file as it may be useful during the picking process.
# 3. Identify the target interval (where the bifurcation is located within the dataset). Set the top and basal locations.
# 4. Identify the maximum and minimum frequencies you would like to search for birfuraction peaks within. The minimum should be lower and the maximum 
#    should be higher than the bifurcation, respectively.
# 5. Check you have the right sampled frequencies and the correct target depth by running a second "high resolution" EHA. Use this EHA to identify target and error frequencies.
# 6. Run the code line by line carefully where instructed and the carefully read the instructions below for the picking of the amplitudes.
# 7. Complete the amplitude difference matrix and calculate your hiatus duration and error duration. This will also give you the error associated in that pick as a percentage.
#
# The base code used for this and the package it is built within is Astrochron for R (Myers, 2014). This code is simply using aspects of Astrochron and streamlining the process of bifurcation picking. Full credit a citation should be as follows:
# Meyers, S.R. (2014). Astrochron: An R Package for Astrochronology. https://cran.r-project.org/package=astrochron
