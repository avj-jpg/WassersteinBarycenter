# Explanation of command-line arguments:
    # -np:      Number of mpi cores (power of 2: 8, 16, 32, etc.).
    # -img:     Text file containing names of input objects.
    # -tC:      Testcase: 1 for Gaussians 2 for everything else.
    # -N:       Number of species (densities) in the system.
    # -nt:      Number of coarse time elements.
    # -rs:      Number of times to refine the mesh uniformly in serial (refines spatial mesh).
    # -rp:      Number of times to refine the mesh uniformly in parallel (refines spatial mesh).
    # -gr:      Number of geometric refinements done prior to order refinements (refines spacetime mesh).
    # -or:      Number of order refinements. Polynomial used are of degree 2^{or}
    # -alpha:   Reaction strength.
    # -beta:    Interaction strength.
    # -alg:     Number of PDHG iterations.
    # -tolPDHG: PDHG tolerance.
    # -cgI:     Number of CG iterations.
    # -ptI:     Number of Brent iterations.
    # -rMax:    Maximum density value for brent solver.
    # -nx0:     Resolution of voxel data (Not used in Testcase 1).
    # -pv:      Save data files for ParaView visualization (-no-pv otherwise).
    # -iP:      Output print interval (-iP 10 prints every 10 steps).

mpirun                                    \
    -np 32                                \
    -bind-to core:2                       \
    ../bary_pdhg                          \
    -img doubleTorusBunny.txt             \
    -tC 2                                 \
    -N 2                                  \
    -nt 1                                 \
    -rs 1                             	  \
    -rp 0                           	  \
    -gr 2                           	  \
    -or 2                           	  \
    -alpha 0                	    	  \
    -beta 0.001                       	  \
    -alg 1000                       	  \
    -tolPDHG 0.000001               	  \
    -cgI 1                           	  \
    -ptI 20                         	  \
    -rMax 40.0                      	  \
    -nx0 256                        	  \
    -iP 50                          	  \
    -pv    
