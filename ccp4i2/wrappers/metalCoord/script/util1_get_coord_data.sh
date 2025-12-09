#!/bin/bash

# Iterate from 1 to 17 with leading zeros
for i in {2..17}; do
    # Leading zeros
    i_zeros=$(printf "%02d" $i)
    echo ${metal} $i ${i_zeros}
    metalCoord coord -n $i --cod -o ${i_zeros}.json
done

metalCoord coord -n 24 --cod -o 24.json
