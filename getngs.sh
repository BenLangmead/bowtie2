#!/bin/sh

NGS_VER=2.9.2

if [ ! -f "ngs-${NGS_VER}/install/include/ngs/Alignment.hpp" ] ; then
    if [ ! -d "ngs-${NGS_VER}/ngs-sdk" ] ; then
        wget https://github.com/ncbi/ngs/archive/${NGS_VER}.tar.gz
        tar zxvf ${NGS_VER}.tar.gz
        rm -f ${NGS_VER}.tar.gz
    fi
    cd ngs-${NGS_VER} && ./configure --prefix=`pwd`/install && make && make install
fi

VDB_VER=2.9.2-1

if [ ! -f "ncbi-vdb-${VDB_VER}/install/include/ncbi-vdb/NGS.hpp" ] ; then
    if [ ! -d "ncbi-vdb-${VDB_VER}/vdb3" ] ; then
        wget https://github.com/ncbi/ncbi-vdb/archive/${VDB_VER}.tar.gz
        tar zxvf ${VDB_VER}.tar.gz
        rm -f ${VDB_VER}.tar.gz
    fi
    cd ncbi-vdb-${VDB_VER} && ./configure --prefix=`pwd`/install --with-ngs-sdk=../ngs-${NGS_VER}/install && make && make install
fi
