#!/bin/sh

d=`dirname $0`
for i in `ls $d/*.pm` ; do echo $i ; perl $i --test ; done
