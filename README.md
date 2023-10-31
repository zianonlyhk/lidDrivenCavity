# Solving the lid-driven cavity problem

## Table of Contents

- [About](#about)
- [Getting Started](#getting_started)
- [Usage](#usage)

## About <a name = "about"></a>

This is a coding practice I did as part of my MPhil at the University of Cambridge. It marks my progress learning C++ and other tools in the developer suite.

## Getting Started <a name = "getting_started"></a>

Simply clone the repo and put it somewhere in your computer. Things are portable.

### Prerequisites

C++ compiler and Make software are required but they are shipped by default in most of the OS out there. The only thing special is `gnuplot` if one wants to use the codes to the fullest extent.

## Usage <a name = "usage"></a>

I first edit the absolute directory of the repo on line 30 of the file "`config`" with the correct value, in this case it's "`/Users/zian/Room214N/github/lidDrivenCavity/`". Notice that it ends with a "`/`". One might also need to make changes to the include directory of eigen in `Makefile`, and the directory to the gnuplot binary in `gnuplot.gp`.

I will then type the following command on my mac terminal under the repo directory:

````
make && ./bin/main && ./gnuplot.gp
````

An example set of x-velocity solutions is presented as:

![](https://raw.githubusercontent.com/zianonlyhk/lidDrivenCavity/dev/uPlot.gif)
