Python script for automating the creation of NACA airfoil

* plot
* export as txt from TE to LE and LE to TE
* export gmsh .geo file for c-mesh generation

Uses generic equations for NACA 4 digits. For 6-series uses the software naca456 from https://www.pdas.com/naca456pdas.html, which is written in fortran and must be compiled before execution of the script.
