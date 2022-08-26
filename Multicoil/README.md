Folder contains codes for TCI submission 2022. (ktslr, L+S-Otazo, L+S-Lin author provided code is modified for reconstruction of single coil (instead of multi coil). We used same parameters as author provided.

=================================================================================

Comparison of ktslr, L+S-Otazo, L+S-Lin, altGDMin-MRI and altGDMin-MRI2.

To generate Table III results:

1. This code requires the Matlab version of the Michigan Image Reconstruction Toolbox (MIRT) from [http://web.eecs.umich.edu/~fessler/code/index.html]
Please set up MIRT before running the examples.
Run the mirt-main/setup.m: L+S-Lin code requires the Matlab version of the Michigan Image Reconstruction Toolbox (MIRT).



Run Main_files_retrospective.m: This run all the 5 algorithms and calculate NMSE (Normalized Mean Square Error), Time required and Similarity index and save these results in Comparison_error.txt [Error(Time)] amd Comparison_sim.txt [sim(Time)] for Brain data.

===================================================================================

 For other datasets and questions contact email sbabu@iastate.edu, namrata@iastate.edu.



If you are using our code please cite our paper: S. Babu, S. S. Nayer, S. G. Lingala and N. Vaswani, "Fast Low Rank Column-Wise Compressive Sensing For Accelerated Dynamic MRI," ICASSP 2022 - 2022 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), 2022, pp. 1346-1350, doi: 10.1109/ICASSP43922.2022.9747549.

