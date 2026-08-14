# User guide

## Setup 


1. Clone the qci library:

    ```bash
    git clone git@github.com:koroeder/QCI-library.git
    ```


2. Load/install libraries

   
   **Ubuntu**

   Install missing libraries with: `sudo apt install liblapack-dev libopenblas-dev gfortran make cmake`

3. Create build directory in the standaloneQCI dir

    ```bash
    cd QCI-library/standaoneQCI      
    mkdir build && cd build                
    ```


4. Run cmake

    ```bash
    FC=gfortran cmake ..
    ```


5. Run make
    ```bash
    make
    ```

## Usage 

* Files needed (Amber): 
  * `start` & `finish` - files containing only xyz coordinates, same atom order
  * `perm.allow` - permutations file (script generated)
  * `coords.prmtop` - Amber topology file
  * `QCI_params.dat` - QCI parameter file
* Running QCI:
  * `./QCI QCI_params.dat > output` 
* Output files:
  * `output` -  logfile
  * `int.xyz`  - interpolated images in trajectory format
  * `int.EofS` - energy band 

## File prep

### `start` type files

* Script  `QCI-library/scripts/inpcrd2start.py` can be used to convert Amber coordinate file to `start` type file
* Use: `python3 inprcd2start.py coords.inpcrd start`

### `perm.allow`

* Script used to create `perm.allow` file is available in the `scripts` dir. 


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
* Plot from output
  * Total energy vs Step 

    ```bash
    awk '/E total/{print $4}' output | gnuplot -p -e "set xlabel 'Step'; set ylabel 'Energy'; plot '-' u 1 w p pt 6; set nokey"
    ```
  * Column Energy type mapping in the output:
    * E total: 4 
    * RMS: 6
    * E rep: 9
    * E constr: 12
  * For E spring:

    ```bash
    awk '/E spring/{print $3}' output | gnuplot -p -e "set xlabel 'Step'; set ylabel 'Energy'; plot '-' u 1 w p pt 6; set nokey"
    ```
    * Use column 6 for E dih