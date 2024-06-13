# Simple Molecular Dynamics Simulation Program

This is a basic molecular dynamics simulation program designed to illustrate how to write a fundamental Molecular Dynamics simulation. It includes features such as the cell method (scaling with O(n log n)), temperature and pressure coupling, many-body interactions, and more.

The code is largely based on the book titled *The Art of Molecular Dynamics Simulation*.

[The Art of Molecular Dynamics Simulation](https://www.cambridge.org/core/books/art-of-molecular-dynamics-simulation/57D40C5ECE9B7EA17C0E77E7754F5874)

## Disclaimer
This is a very old project with many "bad" coding practices!

## Installation
Use the typical CMake installation process.

## How to Run
To start a new simulation, run the following command:
```./dpd```

3 files are provided to start a new simulation
1. dpd.init: define MD parameters, such as temperature, cell size, etc.
2. test.top: define topology of a polymer chain
3. s10000.str: define the starting coordinates.