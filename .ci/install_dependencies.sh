#!/usr/bin/env bash
set -vexu

install_root=$1

apt-get update
apt-get install -y software-properties-common
apt-add-repository universe
apt-get update

apt-get install -y \
  build-essential \
  cmake \
  automake \
  gcc \
  gdb \
  git \
  python3 \
  python3-pip \
  python3-pysam \
  python3-setuptools \
  wget \
  zlib1g-dev

pip3 install tox

if [ ! -d $install_root ]; then
  mkdir $install_root
fi
cd $install_root

#________________________ minimap2 __________________________#
cd $install_root
MINIMAP2_V=2.24
wget https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_V}/minimap2-${MINIMAP2_V}.tar.bz2
tar xf minimap2-${MINIMAP2_V}.tar.bz2
rm minimap2-${MINIMAP2_V}.tar.bz2
cd minimap2-${MINIMAP2_V}
arch_is_arm=$(dpkg --print-architecture | grep '^arm' | wc -l)
if [[ $arch_is_arm -gt 0 ]]
then
    make arm_neon=1 aarch64=1
else
    make
fi
cd ..
cp -s minimap2-${MINIMAP2_V}/minimap2 .

#________________________ mummer ____________________________#
cd $install_root
wget -q https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
tar -xvf mummer-4.0.0rc1.tar.gz
cd mummer-4.0.0rc1
./configure LDFLAGS=-static
make
make install
ldconfig
cd ..

#________________________ racon _____________________________#
git clone --recursive https://github.com/lbcb-sci/racon.git racon-git
cd racon-git
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ../../
cp -s racon-git/build/bin/racon .

