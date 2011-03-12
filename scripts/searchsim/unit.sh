#!/bin/sh

for i in `ls *.pm` ; do echo $i ; perl $i ; done
