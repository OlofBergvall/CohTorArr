# Cohomology-of-moduli-spaces-of-Del-Pezzo-surfaces
SageMath code for computations related to the cohomology of moduli spaces of geometrically marked Del Pezzo surfaces of degrees 3 and 4.

The simplest way to use the program is the following:

1. Load the files into Sage.
2. Initiate a root system, e.g. "ra2=RootSystem(["A",2])".
3. Compute the Weyl group, e.g. "wa2=ra2.ambient_space().weyl_group()".
4. Generate representatives for the conjugacy classes of the Weyl group, e.g. "ca2=wa2.conjugacy_classes_representatives()".
5. Compute the values of the equivariant Poincaré polynomial for all conjugacy classes, e.g. "polsa2=many_pols(ca2,ra2,0,len(ca2)-1)".

There is functionality for computing e.g. the value at a single conjugacy class or computing the Poincaré polynomial for
the action extended by {+/-1}. More information about this is provided in the comments (you might also want to consult the code for the
method "many_pols" in "toric_generate_polynomial.sage").
