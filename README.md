# GooFit used by GooStats
Building status of main: ![main](https://github.com/GooStats/GooFit/actions/workflows/autoTest.yml/badge.svg?branch=main)

GooStats for Borexino

## installation guide

```
git clone --recursive git@github.com:GooStats/GooFit

## then salloc ...
mkdir build_GooFit
cd build_GooFit
cmake ../GooFit
make -j
```
that's it.

## test
```
cd build_GooFit
ctest
```
