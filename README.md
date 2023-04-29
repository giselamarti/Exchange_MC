# Exchange Monte Carlo (Parallel tempering) for an Ising spin glass

The project is divided into 3 folders; the first one for parts 1 and 2, the second one for part 3 and the third one for parts 4 and 5. Each folder contains the 
code and the input and output files for the respective part. The code has been divided into different parts for each part of the project, since different
variables are used, and the general structure of the code is also changed slightly in order to adapt to the requirements of each part.

All the code has been written if fortran90 language and it is executed using the gfortran compiler.

All the plots shown in the written report are computed with Pyhton's library matplotlib, by importing the output files obtained in the main program. 

The file "gaussian_numbers.txt" is used in all the programs and it has been generated using the program gauss.f, with seed=1234, sigma=1, nsamples=1000.
