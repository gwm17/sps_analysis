#sps_analysis
Program to perform offline analysis of Super-Enge Splitpole data at FSU using ROOT. Takes a raw data input .root file from sps evt2root and returns a histogram .root file, a corrected .root file, and a cleaned (no background).

The included Makefile will compile the program. There are three modes in which to run the analysis:
1. -r <dataname>
2. -f <dataname>
3. -a <dataname>
4. -b <dataname>
The dataname refers to the name of the root file without the .root extension

Also contains a program peakfit that fits a range of peaks in the spectrum.

#Usage: 
As the analysis program loops throgh the data it will ask for the user to make cuts/slices on the data. The program will allow the user to manipulate the plot shown for cutting as need until the user double clicks on the canvas. For a 1D cut the user will be prompted to enter max/min values for the histogram on the command line. 

--Example-- 

Enter tsum min: 3450
Enter tsum max: 3700

For a 2D cut, before double clicking use the TCutG tool located in the toolbar to draw a closed shape around the area of interest. After the cut is made double click to move on. Note that double clicking before the cut is complete will result in a crash. After making a succesful cut the canvas will either go black if there are more cuts to be made, or close/go blank if there are no further cuts. If it is black, just click once more on the canvas to bring up the next histogram. 

-r tells the program to run all of the analysis tools. It will loop through the data and make the standard analysis histograms (written to dataname_histo.root) and then performs aberration corrections (written to dataname_corr.root). 
-f only runs the aberration corrections. Note that this requires that the standard histogram file already be created and properly filled.
-a only runs the standard analysis.
-b runs background removal. Requires that the corrected file has already been created and filled.

Standard analysis sorts the data into a cleaned spectrum which contains only the data of interest.
Aberration correction takes takes the x|theta information and corrects away the leading order terms by fitting 3rd order polynomials to well defined peaks in the data and interpolating across the entire set. The correction method will ask the user to input how many polynomials are to be made. A standard number is around 5 polynomials, which should be spread across the entire width of the focal plane detector. 
If the user needs background removal, the Backgnd class will estimate the background using ROOT's TSpectrum tool, and produce a cleaned spectrum. Currently, the program expects the corrected file to be fed for cleaning, and that the corrected position spectrum is the specific spectrum to be cleaned

The peakfit program takes in a ROOT file with histograms and then asks the user to supply ranges to perform a fit over for multiple functions (gaussians, breit-wigner, etc.). It first fits each individual peak and then uses the parameters from the individual fits as a initial guess for the parameters of a global fit. It will then save the results of the fit in a txt file specified by the user.
#Execution:

make

Analysis code:
./analysis -r <dataName> OR
./analysis -f <dataName> OR (only if histogram file has already been created)
./analysis -a <dataName> OR
./analysis -b <dataName>    (only if corrected file has already been created)

Peakfit code:
./peakfit <inputfile> <outputfile>

#Requirements:
ROOT ver. 5 (or newer)
c++11
Raw data stored in .root files from sps evt2root (otherwise tree branch names may need edits)

If you encounter any bugs/see areas of improvement, let me know at gwm17@my.fsu.edu

