# WassersteinBarycenter
This README file provides an overview of the repository and instructions on how to build and use the provided MFEM code to compute the generalized Wasserstein barycenter of multidensity systems. 

# Installation
To use the code in this repository, you need to install and set up the necessary dependencies, with the core requirement being the parallel version of the MFEM library. We will use the jkop-mfem branch of the MFEM fork available [here](https://github.com/pazner/mfem/tree/mfem-jkop). Building the parallel version of MFEM requires an MPI C++ compiler and external libraries such as hypre and METIS. Below are the instructions to download and install MFEM along with its prerequisites. A more detailed version can be found [here](https://mfem.org/building/#parallel-mpi-version-of-mfem).

  1. Create an empty directory 'MFEM/' and set as current working directory:
     ```
     mkdir MFEM
     cd MFEM
     ```
  2. Download and build hypre:
     ```
     git clone https://github.com/hypre-space/hypre.git 
     cd hypre/src
     ./configure && make install
     cd ../../
     ```
  3. Download and build METIS (4.0.3)
     ```
     wget https://github.com/mfem/tpls/raw/gh-pages/metis-4.0.3.tar.gz
     tar -zxvf metis-4.0.3.tar.gz 
     cd metis-4.0.3
     make OPTFLAGS=-Wno-error=implicit-function-declaration
     cd ../
     ln -s metis-4.0.3 metis-4.0
     ```
  4. Download and build the mfem-jkop branch
     ```
     git clone --single-branch --branch mfem-jkop https://github.com/pazner/mfem.git
     cd mfem
     make parallel -j 4 MFEM_USE_ZLIB=YES
     cd ../
     ```
Having successfully built all the dependencies, you can proceed to build the code in this repository by following these instructions:
 1. Clone the repository
    ```
    git clone https://github.com/avj-jpg/WassersteinBarycenter.git
    cd WassersteinBarycenter
    ```
 2. Set the MFEM path in Make.user file:
    ```
    MFEM_DIR = ../mfem
    ```
    (Create the file if it does not exist and set MFEM_DIR)
 3. Compile the code
    ```
    make
    ```
# Usage
From within the `WassersteinBarycenter` directory, you can run a barycenter computation of Gaussian species using:
```
mpirun -np 16 -bind-to core:2 ./bary_pdhg -N 3  -or 0 -gr 1 -rs 1 -alg 100 -cgI 20
```
While a similar barycenter on a 2D embedded surface can be run using:
```
mpirun -np 16 -bind-to core:2 ./surf_pdhg -N 3 -m data/shell_quad.mesh -or 0 -cgI 50 -tC 2 
```
For a detailed explanation of the commandline arguments please refer to one of the example bash scripts found in `examples` or run `bary_pdhg -h`.

# Examples
The numerical examples from our paper can be found in the `examples/` directory. The following commands will run the corresponding test:
 1. Three Gaussians
    ```
    cd examples
    ./gaussians.sh
    ```
 2. Torus and double-torus
    ```
    ./torusDoubleTorus.sh
    ```
 3. Double-torus and bunny
    ```
    ./doubleTorusBunny.sh
    ```
 4. Torus and bunny
    ```
    ./torusBunny.sh
    ```
 5. Torus, double-torus and bunny
    ```
    ./torusDoubleTorusBunny.sh
    ```
 6. Surface Gaussians
    ```
    ./surfaceGaussians.sh
    ```
 7. Indicator functions on a dinosaur surface mesh
    ```
    ./dino.sh
    ```
Reactions and interactions may be added by editing the above script files and changing the `-alpha` and `-beta` flags.

# License
This project is licensed under the MIT License - see the LICENSE file for details.
The MIT License is a permissive license that allows you to use, modify, and distribute the code for both commercial and non-commercial purposes. You are free to do so as long as you include the original copyright notice and disclaimer in all copies or substantial portions of the software.
For more information about the MIT License, please visit the Open Source Initiative or refer to the LICENSE file included with this project.
