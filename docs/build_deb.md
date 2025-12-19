# Build on Debian

0. Install missing packages
```
sudo apt install liblapack-dev libopenblas-dev gfortran make cmake
```
1. Download the source code

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
5. Usage: 
```
./QCI <n_atoms> <params_file>
```

