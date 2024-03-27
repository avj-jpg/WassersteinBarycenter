# Explanation of command-line arguments:
    # -np:      Number of mpi cores (power of 2: 8, 16, 32, etc.).
    # -tC:      Testcase: 1 for Gaussians 2 for everything else.
    # -N:       Number of species (densities) in the system.
    # -m:       Path to Mesh file.
    # -nt:      Number of coarse time elements.
    # -rs:      Number of times to refine the mesh uniformly in serial (refines spatial mesh).
    # -rp:      Number of times to refine the mesh uniformly in parallel (refines spatial mesh).
    # -gr:      Number of geometric refinements done prior to order refinements (refines spacetime mesh).
    # -or:      Number of order refinements. Polynomial used are of degree 2^{or}
    # -s: 	Scale mesh by a factor of 0.01 (-no-s for no scaling).
    # -alpha:   Reaction strength.
    # -beta:    Interaction strength.
    # -alg:     Number of PDHG iterations.
    # -tolPDHG: PDHG tolerance.
    # -cgI:     Number of CG iterations.
    # -ptI:     Number of Brent iterations.
    # -rMax:    Maximum density value for brent solver.
    # -pv:      Save data files for ParaView visualization (-no-pv otherwise).
    # -iP:      Output print interval (-iP 10 prints every 10 steps).


mpirun                              \
    -np 32                          \
    -bind-to core:2                 \
    ../surf_pdhg                    \
    -tC 4                           \
    -N 3                            \
    -m "../data/dino0.vtk"	    \
    -nt 4                           \
    -rs 0                           \
    -rp 0                           \
    -gr 0                           \
    -or 2                           \
    -s				    \
    -alpha 0.1                      \
    -beta 0                         \
    -alg 10000                      \
    -tolPDHG 0.001                  \
    -cgI 100                        \
    -ptI 100                        \
    -rMax 40.0                      \
    -iP 10                          \
    -pv                     
    


