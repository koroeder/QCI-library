# Parameter options


## Admin 

* `COMMENT` 
  * Make a comment
  * Example: `COMMENT this is a comment`
* `DEBUG`
  * Get extended output 
  * Example: `DEBUG`

## General Interpolation settings 

* `MAXITERATIONS <int>`
  * Maximum number of iteration steps. Needs to be at least (preferably greater) than number of atoms. 
  * Units: None  
  * Example: `MAXITERATIONS 100`


* `QCIREADGUESS <filename>`
  * Start from restart/guess 
  * Example: `QCIREADGUESS guess.xyz` 

## Topology options 

`QCIMODE <Topology_type>`

* Topology type  
* Valid types: `AMBER,` `HIRE`, `SBM`, `GEOMETRY`
* Example: `QCIMODE AMBE`

### Amber  

* `TOPFILENAME <filename>`
  * Amber topology file  
  * Example: `TOPFILENAME coords.prmtop`  
* `AMBERCONSTRFILE <filename>`
  * Amber constraints file  
  * Example: 
* `DETECTBASEPAIRS`
  * Automatically detect base-pairs and add base-pair constraints 
  * Example:  `DETECTBASEPAIRS`

### Hire  

* `HIRETOPFILE <filename>`
  * Hire topology file 
  * Example: 
* `HIRECONSTRFILE <filename>`
  * Hire constraints file  
  * Example:

### SBM     

* `SBMCONTACTFILE <filename>`
  * SBM topology file  
  * Example: 

### Geometry 

* `GEOMFILE <filename>` 
  * Geometry file  
  * Example: 

## Atom adding options 

* Options for coordinate calculation algorithms - **must choose one**: 
  * `TRILATERATE` - *(recommended)* use intersection of three spheres to calculate next atom position 
  * `TRILATERATE2` - *(experimental)* extension for `TRILATERATE` for cases where we fail to find three sphere intersection 
  * `USEINTERNALS` - use internal coordinates to calculate next atom position 
  * `USEFOURATOMS` - use for four atom basis for interpolation 
* Strategies for the activation order : 
  * `QCIDOBACK` - Add backbone atoms first 
  * `QCIDOBACKALL` - Add the whole backbone at once 
  * `QCIADDACID`- When residue is started it has to be added entirely. Default for Amber, not compatible with `QCIDOBACK`
* `MINIMISEAFTERADD <int>` 
  * Minimise for N steps after adding next atom. Significantly improves convergence, but slows down interpolation 
  * Units: None  
  * Example: `MINIMISEAFTERADD 10` 
* `DISTCUTADDATOM <double>`
  * Distance cutoff for creating local axis system 
  * Units: Angstrom
  * Example: `DISTCUTADDATOM 50.0` 
* `MAXSEPADDATOM <int>`
  * How far do we look in index separation 
  * Units: None
  * Example: `MAXSEPADDATOM 15`

## Penalty functions 

* `INTMINFACTOR <double>`
  * Scaling factor for internal minima for both repulsion and constraints 
  * Units: None  
  * Example: `INTMINFACTOR 1.0`

### **Repulsion**

* `K_REP <double>` 
  *  Strength of the repulsion term.
  * Units: energy 
  * Example: `K_REP 0.2`
* `REPULSIONCUTOFF <double>`
  * Minimum repulsion cutoff. Used if `DMIN` < `REPULSIONCUTOFF`, where `DMIN` is distance between two atoms for which we calculate the repulsion  
  * Units: Angstrom 
  * Example: `REPULSIONCUTOFF 11.0` 
* `CHECKREPCUTOFF <double>`
  * Factor for checking repulsion list `CHECKREPCUTOFF*REPCUT`
  * Units: None 
  * Example: `CHECKREPCUTOFF 1.25`
* `CHECKREPINTERVAL <int>`
  * Repulsion list update frequency 
  * Units: None 
  * Example: `CHECKREPINTERVAL 10`
* `REPINTMINSEP <int>` 
  * Minimum separation of atom indices to calculate repulsion between them 
  * Units: None 
  * Example: `REPINTMINSEP 4`
  * *Note: Optim default is 20* 

### Constraints 

* `K_CONST <double>`
  * Strength of the constraint term 
  * *Units: Energy* 
  * `K_CONST 1.0D0` 
* `CHECKINTMINCONSTR`
  * Check internal minima for constraints 
  * Logical 
  * Example: `CHECKINTMINCONSTR`
* `MAXCONUSE <int>` 
  * Maximum number of constraints per atom 
  * Units: None 
  * Example: `MAXCONUSE 10` 
* `CONCUTABS <double>`
  * Constraints cutoff ($r_{cutoff}$}). If interpolation gets stuck, we increase this parameter. See 'Potential functions' for more details.  
  * Units: Angstrom 
  * Example: `CONCUTABS 0.1` 
* `CONCUTFRAC <double>`
  * Use fraction for increasing defining $r_{cutoff}$. Alternative for `CONCUTABS` 
  * Units: None
  * Example: `CONCUTFRAC 0.25`
* `CONACTINACT <double>` 
  * Fraction for scaling active-inactive constraints
  * Units: None
  * Example: `CONACTINACT 0.1`

#### **Congeom ONLY settings** 

* `QCICONSTRAINTTOL <double>` 
  * Constraints tolerance used only for creation of geometric constraints 
  * Units: Angstrom? 
  * Example: `QCICONSTRAINTTOL 1.0D-3`


* `QCICONSEP <int>`
  * Max separation in sequence to have geometric constraints 
  * Units: None  
  * Example: `QCICONSEP 15`
* `QCICONCUT <double>`
  * Maximum distance between two atoms to create a constraint 
  * Units: Angstrom 
  * Example: `QCICONCUT 6.0` 

  \

### Dihedrals   

* `DIHEDRALCONSTR <double>` 
  * Strength of the dihedral potential  
  * Units:   
  * `DIHEDRALCONSTR 10.0`
* `DIHEDRALS <int>`
  * Type of dihedral potential. 
  * `1` = chiral centres only *(default)*, `2`= all dihedrals *(experimental)*
  * Example: `DIHEDRALS 1`

## Spring  

* `KSPRING <double>`
  * Strength of the spring potential  
  * Units: Energy / Angstrom^2
  * Example: `KSPRING 5.0`
* `KSPRINGSCALING <double>`
  * Scaling for spring adjustment 
  * Units: None  
  * Example: `KSPRINGSCALING <double>`
* `ADJUSTSPRING <int>`
  * Frequency of spring adjustment/checks  
  * Units: None  
  * Example: `ADJUSTSPRING 10`
* `KSPACINGDEV <double>`
  * Deviation based tolerance for adjusting (raising or lowering) the spring constant  
  * Units: 
  * Example: `KSPACINGDEV 12.0`
* `KMIN <double>`
  * Minimum strength of spring potential  
  * Units: Energy / Angstrom^2
  * Example: `KMIN 0.5`
* `KMAX <double>` 
  * Maximum strength of spring potential 
  * Units: Energy / Angstrom^2
  * Example: `KMIN 7.5`
* `KADJUSTFRAC <double>`
  * Fraction for dynamic spring adjustment 
  * Units: None 
  * Example: 

## Linear groups 

* `USELINGROUPS`
  * Automatically detect groups of atoms which stay fixed during the interpolation. These are added using rigid body quaternion rotation+translation 
  * *Logical* 
  * Example: `USELINGROUPS`

## Image control 

* `NIMAGES <int>`
  * Number of images  
  * Units: None
  * Example: `NIMAGES 10`
* `MAXINTIMAGES <int>`
  * Max number of images that can be added  
  * Units: None  
  * Example: `MAXINTIMAGES 10`
* `MAXIMAGESEPARATION <double>`
  * Image separation that can trigger image addition if there is less than `MAXINTIMAGES`
  * Units: Angstrom / N_atoms  
  * Example: `MAXIMAGESEPARATION 1.0`
* `MINIMAGESEPARATION <double>`
  * Min image separation that triggers image removal  
  * Units: Angstrom / N_atoms  
  * Example: `MINIMAGESEPARATION 0.0`

## Permutations  

* `QCIPERMCHECK <int>`
  * Frequency of the per mutational checks 
  * Units: None  
  * Example: `QCIPERMCHECK 151`
* `QCIPERMCUT <double>`
  * Alignment threshold for permutations. Permutation not applied if alignment is worse in terms of total Euclidean distance  
  * Units: Angstrom  
  * Example: `QCIPERMCUT 1.2` 

## Energy minimisation options  

* `MAXGRADCOMP <double>`
  * Maximum gradient component. Triggers gradient capping. Off by default.
  * Units:  Energy / Angstrom 
  * Example: `MAXGRADCOMP 2.5`


* `DIAGGUESS <double>`
  * Initial guess for diagonal inverse Hessian elements. 
  * Units: None  
  * Example: `DIAGGUESS 0.1`
* `UPDATES <int>`
  * Number of stored information about the previous step  
  * Units: None
  * Example: `DIAGGUESS 0.1`
* `MAXQCIBFGS <double>`
  * Used to take the minimum scale factor for all images for LBFGS step to avoid discontinuities
  * Units: 
  * Example `MAXQCIBFGS 0.2`
* `COLDFUSIONLIMIT <double>`
  * Cold fusion detection criterion
  * Units: Energy
  * Example: `COLDFUSIONLIMIT -1000000000000`

## Convergence criteria 

* `MAXCONSTRAINTE <double>`
  * Max constraint or repulsion energy of any atom
  * Units: Energy
  * Example: `MAXCONSTRAINTE 1.0`
* `RMSTOLERANCE <double>`
  * Max constraint or repulsion force on any atom
  * Units: Energy/Angstrom
  * Example: `RMSTOLERANCE 1.5`
* `MAXERISE <double>`
  * Maximum energy rise from previous step
  * Unites: Energy
  * Example: `MAXERISE 1.0D100`
* `SPRING_GRAD_CONV <double>`
  * Maximum spring force on any atom
  * Units: Energy / Angstrom
  * Example: `SPRING_GRAD_CONV 2.5`
* `QCIRESET <int>`
  *  How often do we check the interpolation progress and adjust convergence criteria
  * Units: None
  * Example:  `QCIRESET 500`

## Output control 

* `DUMPXYZ <int> `
  * Output  frequency
  * Units: None
  * Example: `DUMPXYZ 500`

## Options in development 

* `DETECTBBCROSSINGS <int>` 
  * **does not work atm!** 
  * Frequency for bond crossing detection. Off by default
  * Units: None
  * Example: `DETECTBBCROSSINGS 10`
* `USEIMAGEDENSITY <double>`
  * **does not work atm!** 
  * Image density per unit distance
  * Units: per Angstrom
  * Example: `USEIMAGEDENSITY 10`

## Candidates for removal 

* `QCILINEAR`
* `LINEARCUTOFF`
* `LINEARBB`
* `FREEZEFILE`
* `NMINUNFROZEN`
* `FREEZEFILE`
* `QCIFREEZE`
* `INTCONSTRAINTDEL` - trying to use this keyword terminates QCI
* `QCICONSTRREP` - trying to use this keyword terminates QCI
