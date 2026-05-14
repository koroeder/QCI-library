# Standalone Quasi-Continious Interpolation (QCI) code

### New option to detect backbone crossings
 
 - DETECTBBCROSSINGS int_check_freq

### Note on linear groups - in development
 - Linear groups only work with AMBER atm!
 - Cannot be used together with frozen atoms or QCILINEAR
 - Add parameter option USELINGROUPS to enable


### New paramter options
- K_REP - strength of repulsion potenital (energy units)
- K_CONST - strength of constraint potential (energy units)

### Other important updates:
- Dihedral potential is quadratic and only applied to the chiral centres (DIHEDRALCONSTR).
- Repulsions are removed for any atoms that have a constraint between them (bond, angle, dihedral, etc. ) and for closest 4 atoms in the index sequence. 
- Constraints and repulsions equations were modified from the published ones, so the potential stregth unit is [energy]. Used to be [energy/distance^2].