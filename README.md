# T-Matrix Data Format

This repository contains files for generation and storage of T-matrices in hdf5 format using various software, following the specifications described in a dedicated document.  


 
The example script to retrieve data from a file containing the T-matrix of a TiO2 cylinder is located in the example_data directory together with the hdf5 files.    

The hdf5 files are tracked with git LFS. If these files are not of interest, run

```
git lfs install --skip-smudge 
```

Then clone the repository. 
If they are needed at some point, run
```
git lfs pull
```
If hdf5 files are needed, run
```
git lfs install 
```
and clone the repository.
