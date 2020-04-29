EAMPA Program
#####################################################
Embedded Atom Method Potential Analyser

Written in Python and Fortran (uses F2PY).


What does it do?
#####################################################

Can be used to calculate:
· energy
· energy, forces 
· energy, forces, stress
of a collection of atoms

Calculates some other properties:
· a0
· v0
· e0
· b0
· 9 independent elastic constants (Orthorhombic)
· shear modulus
· young's modulus
· estimated melting temperature



What configuration files can be used?
#####################################################

There is a native format (see examples) similar to an xyz type file.
It will read Quantum Espresso PWscf output files.




What potentials?
#####################################################

EAM potentials.

Works with potentials in the following form:


Pair:

V(r)

Density - one or multiple distinct functions:

rho_1(r), rho_2(r), rho_3(r)

Embedding term - one or multiple terms:

F_1(rho_1) + F_2(rho_2) + F_3(rho_3)


The potentials can be tabulated or analytic.  It's easy to add new, custom functions, either in Python or the F2PY module.






