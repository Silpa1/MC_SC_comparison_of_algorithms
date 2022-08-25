Folder contains codes for TCI submission 2022. 
(ktslr, L+S-Otazo, L+S-Lin author provided code with same parameters is used (cardiac parameters) as author provided.

=================================================================================

Comparison of ktslr, L+S-Otazo, L+S-Lin, altGDMin-MRI and altGDMin-MRI2.

To generate Table III  results:
1. This code requires the Matlab version of the Michigan Image Reconstruction Toolbox (MIRT) from [http://web.eecs.umich.edu/~fessler/code/index.html]
Please set up MIRT before running the examples.
Run the mirt-main/setup.m

2. Load the Datasets (Due to the large size of multicoil dataset cannot be loaded in github. Datasets will be provided on request.) . Copy and paste the .mat files in the same folder contating the files).

3.  Run Main_files_prospective.m: This run all the 5 algorithms and calculate NMSE (Normalized Mean Square Error), Time required and save these results in Comparison_error.txt [Error(Time)].


===================================================================================

For Datasets email sbabu@iastate.edu, namrata@iastate.edu. 

For questions contact sbabu@iastate.edu, namrata@iastate.edu

If you are using our code please cite our paper: S. Babu, S. S. Nayer, S. G. Lingala and N. Vaswani, "Fast Low Rank Column-Wise Compressive Sensing For Accelerated Dynamic MRI," ICASSP 2022 - 2022 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), 2022, pp. 1346-1350, doi: 10.1109/ICASSP43922.2022.9747549.

