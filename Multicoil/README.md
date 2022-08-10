Folder contains codes for TCI submission 2022. 
(ktslr, L+S-Otazo, L+S-Lin author provided code with same parameters is used (cardiac parameters) as author provided.

=================================================================================

Comparison of ktslr, L+S-Otazo, L+S-Lin, altGDMin-MRI and altGDMin-MRI2.

To generate Table III  results:

1.  Run the mirt-main/setup.m: L+S-Lin code requires the Matlab version of the Michigan Image Reconstruction Toolbox (MIRT).

2. Load the Datasets (Due to the large size of multicoil dataset cannot be loaded in github. Datasets will be provided on request.) . Copy and paste the .mat files in the same folder contating the files).

3.  Run Main_files_prospective.m: This run all the 5 algorithms and calculate NMSE (Normalized Mean Square Error), Time required and save these results in Comparison_error.txt [Error(Time)].


===================================================================================

For Datasets email sbabu@iastate.edu, namrata@iastate.edu. 

For questions contact sbabu@iastate.edu, namrata@iastate.edu

If you are using our code please cite our paper: 'Fast Low Rank column-wise Compressive Sensing for Accelerated Dynamic MRI' authr's: Silpa Babu, Sajan Goud Lingala, Namrata Vaswani.

This code is written by Silpa Babu (the code structure is followed from Seyedehsara (Sara) Nayer).
