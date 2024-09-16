# Coupling Matrix Reconfiguration

This project implements the Levenberg-Marquardt optimization algorithm on the orthogonal group for reconfiguring coupling matrices in cross-coupled resonator filter design. It includes several MATLAB functions for complex orthogonalization, random orthogonal matrix generation, and the main optimization routine.

## Files

1. `cm8.mat`: The MATLAB data file storing the canonical coupling matrix `M_can` of an 8-resonator filter and the weight matrix `W` for filter structure.
2. `corth.m`: Performs complex orthogonalization of complex matrices using the Gram-Schmidt method.
3. `rand_corth_mat.m`: Generates a random complex orthogonal matrix.
4. `leven_marq.m`: Implements the Levenberg-Marquardt optimization algorithm on the orthogonal group for coupling matrix reconfiguration.
5. `test_cm_reconfig.m`: A script to test the coupling matrix reconfiguration process.

## Usage

1. Ensure all the `.m` files are in your MATLAB path.
2. Run the `test_cm_reconfig.m` script to see an example of the coupling matrix reconfiguration process.
