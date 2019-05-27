# Usage: ./build [TAG]
ARG TRAVIS_VERSION=packer-1512502276-986baf0
# updated base travis environment: ARG TRAVIS_VERSION=packer-1512502259-986baf0
FROM travisci/ci-garnet:$TRAVIS_VERSION
LABEL Description="Travis CI local build for Zedhat testing"
USER travis
# Install a recent ruby (default is 1.9.3)
RUN bash -lc 'rvm install 2.3.0'
RUN bash -lc 'rvm use 2.3.0 --default'
# Install travis-build to generate a .sh out of .travis.yml
WORKDIR /home/travis/builds
RUN git clone --depth=1 https://github.com/travis-ci/travis-build.git
WORKDIR travis-build
RUN bash -lc 'gem install travis'
# to create ~/.travis
#RUN bash -lc echo "y" | travis
RUN bash -lc 'travis version'
RUN ln -s $(pwd) /home/travis/.travis/travis-build
RUN bash -lc 'bundle install'
# workaround https://github.com/travis-ci/travis-ci/issues/8098
WORKDIR /home/travis/.travis/travis-build
RUN bash -lc 'bundler add travis'
RUN bash -lc 'bundler binstubs travis'
#  Create project dir, assuming your project is `AUTHOR/PROJECT` on GitHub
RUN mkdir /home/travis/builds/boyle
WORKDIR /home/travis/builds/boyle
RUN git clone --depth=1 https://github.com/boyle/zedhat.git
WORKDIR zedhat
# Fix up the missing "sudo" config (apt-get install)
RUN sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
## various 'apt-get update' failures due to outdated GPG keys and dead server certificates
RUN sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 960B2B2623A0BD5D 6B05F25D762E3157 F6609E60DC62814E 0C49F3730359A14518585931BC711F9BA15703C6
#RUN sudo apt-get install --reinstall ca-certificates
#RUN sudo mkdir /usr/local/share/ca-certificates/cacert.org
#RUN sudo wget -P /usr/local/share/ca-certificates/cacert.org http://www.cacert.org/certs/root.crt http://www.cacert.org/certs/class3.crt
#RUN sudo update-ca-certificates 
#RUN grep ^ /etc/apt/sources.list /etc/apt/sources.list.d/*
RUN sudo rm /etc/apt/sources.list.d/hhvm-repository.list
#RUN sudo add-apt-repository --remove ppa:PPA_Name/ppa
RUN sudo apt-get update -q
RUN sudo apt-get  install  gcc-6 gfortran libblas-dev libopenblas-dev libopenblas-base libsuitesparse-dev liblapacke liblapacke-dev lcov tree
# change to the branch or commit you want to investigate
# compile travis script into bash script
RUN bash -lc '~/.travis/travis-build/bin/travis compile > ci.sh'
RUN sed -i 's|\(--branch\\=\)[^ ]* |\1master\\ |' ci.sh
RUN sed -i 's|\(travis_cmd export\\ CC\\=\).*|\1\\"gcc-6\\" --echo|' ci.sh
RUN chmod +x ci.sh
#RUN bash ci.sh

ENTRYPOINT /bin/bash

