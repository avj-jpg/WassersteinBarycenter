# Explanation of command-line arguments:
    # -np:      Number of mpi cores.
    # -nt:      Number of coarse time elements.
    # -rs:      Number of times to refine the mesh uniformly in serial (refines spatial mesh).
    # -rp:      Number of times to refine the mesh uniformly in parallel (refines spatial mesh).
    # -alg:     Number of PDHG iterations.
    # -tolPDHG: PDHG tolerance.
    # -gr:      Number of geometric refinements done prior to order refinements (refines spacetime mesh).
    # -or:      Number of order refinements. Polynomial used are of degree 2^{or}.
    # -N:       Number of species (densities) in the system.
    # -tC:      Testcase. 1 for gaussian and 2 for everything else. 
    # -pv:      (boolean) Save data files for ParaView visualization.




mpirun \
    -np 16 \
    -bind-to core:2 \
    ./bary_pdhg \
    -nt 4 \
    -rs 1 \
    -rp 0 \
    -alg 10000 \
    -tolPDHG 0.000001 \
    -gr 0 \
    -or 0 \
    -N 3  \
    -tC 1 \
    -pv  \
    -iP 1 \
    -alpha 0 \
    -beta 0 \
    -img gaussians.txt 

