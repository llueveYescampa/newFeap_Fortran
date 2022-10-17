#!/bin/bash

for i in $@; do
    echo $i
    cd $i
    ls 
    rename 's/.for$/.f90/' *.for
    ls 
    cd ..
done

