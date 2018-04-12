# zedhat ẑ  [![Build Status](https://travis-ci.org/boyle/zedhat.svg?branch=master)](https://travis-ci.org/boyle/zedhat) [![Coverity Scan Status](https://scan.coverity.com/projects/15229/badge.svg)](https://scan.coverity.com/projects/boyle-zedhat) [![Coverage Status](https://coveralls.io/repos/github/boyle/zedhat/badge.svg?branch=master)](https://coveralls.io/github/boyle/zedhat?branch=master)
A tool for Electrical Impedance Tomography (EIT) and Electrical Resistivity Tomography (ERT).
Current is applied to electrodes on the boundary and voltages are measured at other electrodes.
From the measured voltages, we reconstruct an estimate of the conductivity ẑ over the interior.

## Dependencies

For Ubuntu 16.04, the following packages need to be installed to compile zedhat:

 - pkg-config
 - libopenblas-dev
 - liblapacke-dev
 - libsuitesparse-dev
 - libmatio-dev

and for developers

 - lcov
 - astyle

```sudo apt install pkg-config libopenblas-dev liblapacke-dev libsuitesparse-dev libmatio-dev```

```sudo apt install lcov astyle```
