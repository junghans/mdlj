# Simple (micro)canonical Molecular Dynamics simulation of a Lennard-Jones fluid

This code was mainly written for educational purposes and to test the performance of different neighbor search algorithms

## Algorithms

* Periodic boundary conditions
* 2nd order Verlet integrator
* 1st order Verlet with thermostat
  * Langevin thermostat
  * DPD thermostat
* Neighbor search
  * N square
  * Verlet lists on N square
  * Linked-cell lists
  * Verlet lists on linked-cell lists

## [Documentation](http://junghans.github.io/mdlj/)

## Issues

Report bugs on the [github issues site](https://github.com/junghans/mdlj/issues)


[![Build Status](https://travis-ci.org/junghans/mdlj.svg?branch=master)](https://travis-ci.org/junghans/mdlj)
[![codecov.io](https://codecov.io/github/junghans/mdlj/coverage.svg?branch=master)](https://codecov.io/github/junghans/mdlj?branch=master)
