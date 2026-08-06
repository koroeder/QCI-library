# User guide

## Setup 


1. Clone the qci library:

```bash
git clone git@github.com:koroeder/QCI-library.git
```


2. Load/install libraries

   
   1. **Create**

      Request interactive session for the compute node

```bash
  srun -p roeder_cpu --time 60 --pty /bin/bash -l
```

       Load modules

```bash
  module load gcc/11.4.0-gcc-11.4.0 cmake/3.27.7-gcc-11.4.0 openblas/0.3.24-gcc-11.4.0
```

     b. **Ubuntu**

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
  * `./QCI <N_atoms> QCI_params.dat > output` 
* Output files:
  * `output` -  logfile
  * `int.xyz`  - interpolated images in trajectory format
  * `int.EofS` - energy band 

## File prep

### `start` type files

* Script  `QCI-library/scripts/inpcrd2start.py` can be used to convert Amber coordinate file to `start` type file
* Use: `python3 inprcd2start.py coords.inpcrd start`

### perm.allow

* Script used to create perm.allow files is available with OPTIM code: <https://www-wales.ch.cam.ac.uk/OPTIM/>
  * relative path in OPTIM source code:  `wales/SCRIPTS/make_perm.allow`
  * Use: `python2 perm-pdb.py conf.pdb AMBER`

### Symmetrise  topology file:

* Script is available with OPTIM code: <https://www-wales.ch.cam.ac.uk/OPTIM/>
  * Relative path: `wales/SCRIPTS/AMBER/symmetrise_prmtop`
  * Use: `python2 perm-prmtop.ff19.py coords.prmtop.old coords.prmtop wales/AMBERTOOLS/dat/leap/lib` 
    * *Note: change path to Amber to absolute path* 


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