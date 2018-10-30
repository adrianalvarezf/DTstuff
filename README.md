How to use the cosmic analysis programs
==============
STEP 1: Histograms

     root DTNTuple_runXXXX.root
     root [0]DTTree->Process("histograms_HH_and_HL.C") 
     
This creates a root file (e.g. runXXXX_histograms.root) .

STEP 2: HL/HH

Compile histodiv.C , then use the output of step1:

      ./histodiv.exe 1 1 runXXXX_histograms.root
      
In the first argument 1 means 1 input file, for the second argument 1 means printgifs=true.

This creates hist_div_runXXXX.root , a txt file with the fit values and the gifs (optional).

STEP 3: Compare peak values of two runs

Compile diff.C , then:

      ./diff.exe hist_div_runXXXX.txt hist_div_runYYYY.txt
Make sure to cut out the first two lines of text in the txts before doing this!

How to use the collision analysis program
==============
Obtain t0s with:

     root DTNTuple_runXXXX.root
     root [0]DTTree->Process("t0_fitter.C") 
This creates the txt file with the t0s (fit_t0_runXXXX.txt). Make sure to cut out the first two lines of this txt before you use it for the corrections!

How to use the corrections program
==============
To calculate the corrections, first compile corrections_calculator.C , then:

      ./corrections_calculator.exe hist_div_runXXXX.root fit_t0_runXXXX.txt
      
How to use the efficiencies program
==============
Compile MyEffWithDigis_all.C , then:

      ./MyEffWithDigis_all.exe directorypath  DTNTuple_runXXXX.root
      
Note
==============
To compile anything you can use compilaC like this:
       ./compilaC  MyEffWithDigis_all
      
      
      
      
      
