# dcs_bfi_fitting
Code for fitting BFI to measurements


## steps to install and run
 
Clone or download the repo:
```bash
git clone https://github.com/ScottJordan/dcs_bfi_fitting.git
```

Change directory to be in the repo: ``cd dcs_bfi_fitting``


run the fitting procedure 

```bash
julia --project=@. bfi_fit.jl -d <data-file>.csv -m <mudata-file>.csv -o <output-file>.csv
```

Replace ``<data-file>.csv`` with the name of the file that contains g2 data. This file should have a header with names ``Time,Detector 1 Correlation,Detector 2 Correlation,Detector 3 Correlation,Detector 4 Correlation`` in the first five columns. It also assumes that every 128 rows are seperate data points and can be sperated by a line starting with CPS.

Replace ``<mudata-file>.csv`` with the name of the file that contains the muA and muS data. This file should have the header ``mA,mS`` and data should start on the 3rd row. Each row in this file corresponds on to one section of the 128 data points in the data-file. 

Replace ``<output-file>.csv`` with the desired name of the output file that will contain the sum of squared errors, beta, and BFI from each fit. 