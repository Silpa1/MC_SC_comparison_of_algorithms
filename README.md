Folder contains codes for TCI submission 2022. 
(ktslr, L+S-Otazo, L+S-Lin author provided code is modified for reconstruction of single coil (instead of multi coil). We used same parameters as author provided.

=================================================================================

Comparison of ktslr, L+S-Otazo, L+S-Lin, altGDMin-MRI and altGDMin-MRI2.

To generate Table III and Table V results:

1.  Run the mirt-main/setup.m: L+S-Lin code requires the Matlab version of the Michigan Image Reconstruction Toolbox (MIRT).

2.  Run Main_comparison_of_algorithm.m: This run all the 5 algorithms and calculate NMSE (Normalized Mean Square Error), Time required and Similarity index and save these results in Comparison_error.txt [Error(Time)] amd Comparison_sim.txt [sim(Time)].


===================================================================================
The folder "Dataset" contains the MRI used in this paper (except long-speech sequence because of its larger size). 

