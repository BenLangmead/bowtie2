#!/bin/bash

set -e

while getopts "sb:" opt; do
    case $opt in
        s) use_sra=1 ;;
        b) branch="$OPTARG" ;;
        *) echo "Usage: $0 [-s] [-b <branch_name>]" && exit 1
    esac
done
shift $(($OPTIND - 1))

if [ "$branch" == "" ] ; then
    branch="master"
fi

set -x
yum install -y git zip unzip pandoc

git clone --recursive https://github.com/BenLangmead/bowtie2.git
if [ $? -ne 0 ] ; then
    echo "Unable to clone bowtie2 repo"
    exit 1
fi

cd bowtie2

git branch -a | grep "$branch" 2>&1 > /dev/null
if [ $? -ne 0 ] ; then
    echo "branch '$branch' does not exist"
    exit 1
else
     git checkout "$branch"
fi

if [ $use_sra -eq 1 ] ; then
    # this variant is needed to compile ncbi-vdb
    source /hbb/activate
    yum install -y java-1.8.0-openjdk-devel.$(uname -m)
    make sra-deps
fi

# this variant creates static binaries with PIC
source /hbb_exe_gc_hardened/activate

mkdir /mybin
echo  'res=`echo $@ | sed "s/-L.*$//"`; /opt/rh/devtoolset-8/root/usr/bin/ar $res;' > /mybin/ar
chmod +x /mybin/ar && export PATH=/mybin:$PATH

# this version of pandoc does not support 'smart' filters
sed -i 's/-smart//' Makefile

make static-libs
if [ $? -ne 0 ] ; then
    echo "Unable to build tbb and/or zlib static dependencies"
    exit 1
fi

make -j4 bowtie2-bin-pkg STATIC_BUILD=1 WITH_ZSTD=1 USE_SRA=$use_sra
if [ $? -ne 0 ] ; then
    echo "Unable to create bowtie2 package"
    exit 1
fi
echo "Running libcheck..."
libcheck bowtie2-{align,build,inspect}-*

echo "Running hardening check..."
hardening-check -b bowtie2-{align,build,inspect}-*

echo "Copying binary package"
cp /bowtie2/*.zip /io
