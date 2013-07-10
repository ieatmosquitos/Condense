#!/bin/bash
philes=`ls $1/*.g2o`
for f in ${philes}
  do
    ./condense $f -sl 10 -ml 1;
    ./condense $f -sl 10 -ml 2;
    ./condense $f -sl 10 -ml 5;
    ./condense $f -sl 10 -ml 10;
    ./condense $f -sl 10 -ml 20;
    ./condense $f -sl 10 -ml 50;
    ./condense $f -sl 10 -ml 1000;

    ./condense $f -sl 30 -ml 1;
    ./condense $f -sl 30 -ml 2;
    ./condense $f -sl 30 -ml 5;
    ./condense $f -sl 30 -ml 10;
    ./condense $f -sl 30 -ml 20;
    ./condense $f -sl 30 -ml 50;
    ./condense $f -sl 30 -ml 1000;
    
    ./condense $f -sl 50 -ml 1;
    ./condense $f -sl 50 -ml 2;
    ./condense $f -sl 50 -ml 5;
    ./condense $f -sl 50 -ml 10;
    ./condense $f -sl 50 -ml 20;
    ./condense $f -sl 50 -ml 50;
    ./condense $f -sl 50 -ml 1000;
    
    ./condense $f -sl 70 -ml 1;
    ./condense $f -sl 70 -ml 2;
    ./condense $f -sl 70 -ml 5;
    ./condense $f -sl 70 -ml 10;
    ./condense $f -sl 70 -ml 20;
    ./condense $f -sl 70 -ml 50;
    ./condense $f -sl 70 -ml 1000;
    
    ./condense $f -sl 100 -ml 1;
    ./condense $f -sl 100 -ml 2;
    ./condense $f -sl 100 -ml 5;
    ./condense $f -sl 100 -ml 10;
    ./condense $f -sl 100 -ml 20;
    ./condense $f -sl 100 -ml 50;
    ./condense $f -sl 100 -ml 1000;
  done