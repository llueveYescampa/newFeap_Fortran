#!/bin/bash

for i in $@; do
    
    echo $i, ${i,,}
     mv $i  ${i,,} 
done

