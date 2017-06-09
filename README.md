# zedhat ẑ  [![Build Status](https://travis-ci.com/boyle/zedhat.svg?token=ppJxXrcSY7ZkfdM9JuFM&branch=master)](https://travis-ci.com/boyle/zedhat)
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

```sudo apt install pkg-config libopenblas-dev liblapacke-dev libsuitesparse-dev libmatio-dev```
