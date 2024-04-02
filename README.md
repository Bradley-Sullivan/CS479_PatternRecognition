## CS 479 - Pattern Recognition - Programming Assignment 2

Author : Bradley Sullivan

Date : April 1st, 2024

## Compilation Instructions

Navigate to the project directory and execute

    make

The executable will be created named 'PA2'. To run execute

    ./PA2

## Dependencies

### meschach

This project makes use of the Meschach matrix math library.
To build, navigate to `lib/mesch12b/` and execute `make all`.

A compiled version of meschach must exist and be placed under `lib/` with the name `libmeschach.a`. 

### gnuplot

For graph plotting, an installation of 'gnuplot' is required (version 6.0 or greater). 

This project also makes use of an ANSI C interfacing library `gnuplot_i`. To compile, navigate to `lib/gnuplot_i-master/` and execute `make`. 

The object file `gnuplot_i.o` must be placed under `lib/`

Matrix library: https://homepage.divms.uiowa.edu/~dstewart/meschach/

Plotting library interface: https://github.com/mithodin/gnuplot_i

gnuplot module: http://gnuplot.info/

