Thesis Scripts
=============

Apart from thesis scripts, this repo also contains `compute_stress_mop`,
a LAMMPS compute that can be used to compute the stress on a plane 
using the [Method of Planes](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.52.1627). 

### Usage

    compute _computeName_ _group_ stress/mop _z_ _x_

Here, the plane on which the 
stress is to be computed is assumed to be 
parallel to the xy plane. $z$ indicates
the $z$ coordinate of the plane and $x$
refers to the bin size along the x axis.
The flow is assumed to be invariant along
the y axis and hence each bin on the plane
extends to the size of the box along the y 
axis.


