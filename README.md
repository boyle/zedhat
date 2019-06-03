# zedhat ẑ  [![Build Status](https://travis-ci.org/boyle/zedhat.svg?branch=master)](https://travis-ci.org/boyle/zedhat) [![Coverity Scan Status](https://scan.coverity.com/projects/15229/badge.svg)](https://scan.coverity.com/projects/boyle-zedhat) [![Coverage Status](https://coveralls.io/repos/github/boyle/zedhat/badge.svg?branch=master)](https://coveralls.io/github/boyle/zedhat?branch=master)
A tool for Electrical Impedance Tomography (EIT) and Electrical Resistivity Tomography (ERT).
Current is applied to electrodes on the boundary and voltages are measured at other electrodes.
From the measured voltages, we reconstruct an estimate of the conductivity ẑ over the interior.

## Dependencies

For Ubuntu 16.04, the following packages need to be installed to compile zedhat:

 - pkg-config
 - gfortran
 - libblas-dev
 - libopenblas-dev
 - libopenblas-base
 - libsuitesparse-dev
 - liblapacke
 - liblapacke-dev

and for developers

 - lcov
 - astyle

For macOS 10.14.4, install Xcode from the Mac store, launch Xcode to complete the installation, then from a terminal install [homebrew](https://brew.sh) and
```
brew install autoconf automake pkg-config lapack openblas suite-sparse
export LDFLAGS="-L/usr/local/lib -L/usr/local/opt/openblas/lib"
export CFLAGS="-I/usr/local/include -I/usr/local/opt/openblas/include"
export PKG_CONFIG_PATH="/usr/local/opt/openblas/lib/pkgconfig"
./autogen && ./configure && make
```
(Homebrew `openblas` is necessary since Xcode does not provide lapacke, only
clapack, the older C interface to LAPACK.
Homebrew `suite-sparse` provides UMFPACK and CHOLMOD.)

 For Windows 10, install [chocolatey](https://chocolatey.org/install) as admin, then from regular PowerShell
 ```
 choco -y install mingw
 choco -y install git.install
 ```
 and finally start "Git Bash"
