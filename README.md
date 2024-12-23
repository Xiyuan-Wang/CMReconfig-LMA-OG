# Coupling Matrix Reconfiguration

This project implements the Levenberg-Marquardt optimization algorithm on the orthogonal group for reconfiguring coupling matrices in filter design. The general isospectral flow algorithm is also implemented for algorithm comparison. It includes several MATLAB functions for complex orthogonalization, random orthogonal matrix generation, and testing scripts for the reconfiguration algorithms.

## Files

1. `cm/cm4.mat`, `cm/cm10.mat`, `cm/cm8.mat`: These MATLAB data files contain the canonical coupling matrix `M0` for 4-, 10-, and 8-resonator filters, respectively. Each file also includes a corresponding binary matrix `W`, where nonzero elements indicate the couplings to be reduced.
2. `corth.m`: Performs complex/real orthogonalization of complex/real matrices using the Gram-Schmidt method.
3. `rand_corth_mat.m`: Generates a random complex orthogonal matrix.
4. `leven_marq.m`: Implements the Levenberg-Marquardt optimization algorithm on the orthogonal group for coupling matrix reconfiguration.
5. `gen_iso_flow.m`: Implements the general isospectral flow algorithm for coupling matrix reconfiguration.
6. `test_cm_reconfig.m`: Tests the coupling matrix reconfiguration algorithms.
7. `reconfig_convergence.m`: Tests the performance of the algorithms in terms of convergence.
8. `reconfig_success_count.m`: Tests the performance of the algorithms in terms of successful reconfigurations.
9. `data`: The folder where the results of the testing code are stored.

## Usage

1. Ensure all the `.m` files are in your MATLAB path.
2. Run the `test_cm_reconfig.m` script to see an example of the coupling matrix reconfiguration process.
3. Run the `reconfig_convergence.m` script to see the convergence performance of the reconfiguration algorithms.
4. Run the `reconfig_success_count.m` script to see the results of the successful reconfigurations.
