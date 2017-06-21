#!/bin/bash
CURRENT_DIR=$(dirname $(readlink -f $0))
cd
cd $CURRENT_DIR
tar cvzf $CURRENT_DIR/gcn.mops.tar.gz gcn.mops
R CMD INSTALL $CURRENT_DIR/gcn.mops.tar.gz
rm gcn.mops.tar.gz

DATASET=bamDataRanges21.rds
while true
do
  read -p "Do you want to run cn.mops & gcn.mops with a small sample (~10-15 minutes)? (y/n): " input

  case $input in
   [yY]* ) R -e  "bamDataRanges <- readRDS(\"$CURRENT_DIR/R_scripts/data/$DATASET\"); source(\"$CURRENT_DIR/R_scripts/test_script.R\")"
           break;;

   [nN]* ) exit;;

   * )     echo "Invalid input";;
  esac
done
