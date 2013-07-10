#!/bin/bash
philes=`ls $1/*.g2o`
for f in ${philes}
  do
    ./cond3nse $f -sl 10 -ml 1;
    ./cond3nse $f -sl 10 -ml 2;
    ./cond3nse $f -sl 10 -ml 5;
    ./cond3nse $f -sl 10 -ml 10;
    ./cond3nse $f -sl 10 -ml 20;
    ./cond3nse $f -sl 10 -ml 50;
    ./cond3nse $f -sl 10 -ml 1000;

    ./cond3nse $f -sl 30 -ml 1;
    ./cond3nse $f -sl 30 -ml 2;
    ./cond3nse $f -sl 30 -ml 5;
    ./cond3nse $f -sl 30 -ml 10;
    ./cond3nse $f -sl 30 -ml 20;
    ./cond3nse $f -sl 30 -ml 50;
    ./cond3nse $f -sl 30 -ml 1000;
    
    ./cond3nse $f -sl 50 -ml 1;
    ./cond3nse $f -sl 50 -ml 2;
    ./cond3nse $f -sl 50 -ml 5;
    ./cond3nse $f -sl 50 -ml 10;
    ./cond3nse $f -sl 50 -ml 20;
    ./cond3nse $f -sl 50 -ml 50;
    ./cond3nse $f -sl 50 -ml 1000;
    
    ./cond3nse $f -sl 70 -ml 1;
    ./cond3nse $f -sl 70 -ml 2;
    ./cond3nse $f -sl 70 -ml 5;
    ./cond3nse $f -sl 70 -ml 10;
    ./cond3nse $f -sl 70 -ml 20;
    ./cond3nse $f -sl 70 -ml 50;
    ./cond3nse $f -sl 70 -ml 1000;
    
    ./cond3nse $f -sl 100 -ml 1;
    ./cond3nse $f -sl 100 -ml 2;
    ./cond3nse $f -sl 100 -ml 5;
    ./cond3nse $f -sl 100 -ml 10;
    ./cond3nse $f -sl 100 -ml 20;
    ./cond3nse $f -sl 100 -ml 50;
    ./cond3nse $f -sl 100 -ml 1000;
  done