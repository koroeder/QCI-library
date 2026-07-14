# Build on Debian

0. Install missing packages
    ```
    sudo apt install liblapack-dev libopenblas-dev gfortran make cmake
    ```
1. Download the source code
    ```
    git clone git@github.com:koroeder/QCI-library.git
    ```
2. Create build directory in the standaloneQCI dir
    ```
    cd QCI-library/standaoneQCI    
    mkdir build
    cd build
    ```
3. Run cmake
    ```
    FC=gfortran cmake ..
    ```
4. Run make
    ```
    make
    ```

## Usage 

* Files needed: 
  * `start` & `finish` - file containing only xyz coordinates, same atom order
  * `perm.allow` - permutations file (script generated)
  * `coords.prmtop` - Amber topology file
  * `QCI_params.dat` - QCI parameter file
* Running QCI:
  * `./QCI <N_atoms> <params_file> > output` 
* Output files:
  * `output` -  logfile
  * `int.xyz`  - interpolated images in trajectory format
  * `int.EofS` - energy band 

## Visualising output 

* Image distance and energy profiles can be easily plotted with provided gnuplot scripts (default output) 
  * ```bash
    gnuplot -p plot_image_dist.txt
    ```
  * ```bash
    gnuplot -p plot_energy_profile.txt
    ```
* VMD visualisation
  * ```bash
    vmd -parm7 coords.prmtop -xyz int.xyz
    ```