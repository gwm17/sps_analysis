#sps_analysis
Program to perform offline analysis of Super-Enge Splitpole data at FSU using ROOT. Takes a raw data input .root file from sps evt2root and returns a histogram .root file and potentially a corrected .root file.

The included Makefile will compile the program. There are three modes in which to run the analysis:
1. -r dataname
2. -f dataname
3. -a dataname
The dataname refers to the name of the root file without the .root extension

#Usage: 
As the program loops throgh the data it will ask for the user to make cuts/slices on the data. The program will allow the user to manipulate the plot shown for cutting as need until the user double clicks on the canvas. For a 1D cut the user will be prompted to enter max/min values for the histogram on the command line. 

--Example-- 

Enter tsum min: 3450
Enter tsum max: 3700

For a 2D cut, before double clicking use the TCutG tool located in the toolbar to draw a closed shape around the area of interest. After the cut is made double click to move on. Note that double clicking before the cut is complete will result in a crash. After making a succesful cut the canvas will either go black if there are more cuts to be made, or close/go blank if there are no further cuts. If it is black, just click once more on the canvas to bring up the next histogram. 

-r tells the program to run all of the analysis tools. It will loop through the data and make the standard analysis histograms (written to dataname_histo.root) and then performs aberration corrections (written to dataname_corr.root). 
-f only runs the aberration corrections. Note that this requires that the standard histogram file already be created and properly filled.
-a only runs the standard analysis.

Standard analysis runs through sum time gating, particle ID in E-dE and dE-position, front/rear gating, Si-coincidence in the scattering chamber, and creates a wide number of histograms.
Aberration correction takes takes the x|theta information and corrects away the leading order terms by fitting 3rd order polynomials to well defined peaks in the data and interpolating across the entire set. The correction method will ask the user to input how many polynomials are to be made. A standard number is around 5 polynomials, which should be spread across the entire width of the focal plane detector. 

#Execution:
./analysis -r dataName OR
./analysis -f dataName OR
./analysis -a dataName

#Requirements:
ROOT ver. 5 (or newer)
Raw data stored in .root files from sps evt2root (otherwise tree branch names may need edits)

If you encounter any bugs/see areas of improvement, let me know at gwm17@my.fsu.edu

#Note: 
The output files will be smaller than the (presumably) large data file, but stil could (especially the corrected file) be up to several GB in size. Keep an eye on storage space.
