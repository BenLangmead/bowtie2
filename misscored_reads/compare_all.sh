#!/bin/bash

echo "======================"
echo "===   SENSITIVE   ==="
echo "======================"

echo "----- E2E Low -----"
python compare.py bt2 low archive-master def2_s $1

echo
echo "----- E2E High -----"
python compare.py bt2 high archive-master def2_s $1

echo
echo "----- Local Low -----"
python compare.py bt2loc low archive-master def0_sl $1

echo
echo "----- Local High -----"
python compare.py bt2loc high archive-master def0_sl $1

echo "======================"
echo "=== VERY SENSITIVE ==="
echo "======================"

echo "----- E2E Low -----"
python compare.py bt2 low archive-master def3_vs $1

echo
echo "----- E2E High -----"
python compare.py bt2 high archive-master def3_vs $1

echo
echo "----- Local Low -----"
python compare.py bt2loc low archive-master def0_vsl $1

echo
echo "----- Local High -----"
python compare.py bt2loc high archive-master def0_vsl $1
