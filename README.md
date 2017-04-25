# DTstuff

How to use
==============
STEP 1: Histograms
     root DTNTuple_runXXXX.root
     root [0]DTTree->Process("histograms_HH_and_HL.C") 
This creates a root file e.g. run291222_3600V_histograms.root .

STEP 2: HL/HH
Compile histodiv.C , then
      ./histodiv.exe 1 1 run291222_3600V_histograms.root
(First 1 means 1 file, second 1 means print gifs = true).
This creates hist_div_run291222_3600V.root , a txt file with the fit values and the gifs (optional).
