# Standalone Quasi-Continious Interpolation (QCI) code

See [Quick start](docs/Qucik_start.md) 

For recent changes see [changelog](changelog.md)

Examples are available in the [testing directory](testing/).

## OpenMP version 
This version has not been extensively tested. 
Run cmake with:
``` FC=gfortran cmake -DWITH_OPENMP=ON .. ```

Run:
``` OMP_NUM_THREADS=<n_cores> ./QCI QCI_parms.dat > output ``` 

https://github.com/user-attachments/assets/b5e576a8-1d54-4662-b727-6894592792a0




