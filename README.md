#sps_analysis
Program to perform offline analysis of Super-Enge Splitpole data at FSU using ROOT. Takes a raw data input .root file from sps evt2root and returns a histogram .root file and potentially a corrected .root file.

Use Makefile to build the program. Program takes a minimum of two and a maximum of three file names at execution. Input order is as follows:
1. name of data file 
2. name of new histogram file
3. name of new corrected data file (optional)

#Usage: 
As the program loops throgh the data it will ask for the user to make cuts/slices on the data. The program will allow the user to manipulate the plot shown for cutting as need until the user double clicks on the canvas. For a 1D cut the user will be prompted to enter max/min values for the histogram on the command line. 
--Example-- 
Enter tsum min: 3450
Enter tsum max: 3700
For a 2D cut, before double clicking use the TCutG tool located in the toolbar to draw a closed shape around the area of interest. After the cut is made double click to move on. Note that double clicking before the cut is complete will result in a crash. After making a succesful cut the canvas will either go black if there are more cuts to be made, or close/go blank if there are no further cuts. If it is black, just click once more on the canvas to bring up the next histogram. 
If you give the program a third name at execution it will perform polynomial interpolation corrections for the aberrations in x|theta. Two names will only perform the standard sorting and analysis of data.

#Execution:
1. make clean
2. make
3. ./analysis dataName histoName correctedName

#Requirements:
ROOT ver. 5 (or newer)
Raw data stored in .root files from sps evt2root (otherwise tree branch names may need edits)

If you encounter any bugs/see areas of improvement, let me know at gwm17@my.fsu.edu

#Note: The output files will be smaller than the (presumably) large data file, but stil could (especially the corrected file) be up to several GB in size. Keep an eye on storage space.
