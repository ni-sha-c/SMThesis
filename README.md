Thesis Scripts
=============

Apart from thesis scripts, this repo also contains `compute_stress_mop`,
a LAMMPS compute that can be used to compute the stress on a plane 
using the [Method of Planes](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.52.1627). 

### Installation

To use this compute, add the files `compute_stress_mop.cpp` and 
`compute_stress_mop.h` (available in the `compute_stress_mop` directory) 
 to your LAMMPS `src` directory. Recompile LAMMPS (if `Makefile.list` 
needs to be explicitly edited, add the compute files and then recompile
LAMMPS) 

### Usage

    compute computeName group stress/mop z x

Here, the plane on which the 
stress is to be computed is assumed to be 
parallel to the xy plane. `z` indicates
the z coordinate of the plane and `x`
refers to the bin size along the x axis.
The flow is assumed to be invariant along
the y axis and hence each bin on the plane
extends to the size of the box along the y 
axis.


