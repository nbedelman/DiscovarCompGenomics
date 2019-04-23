#!/bin/bash

miyagiOut=$1
converted=$2

perl -pe 's/],/],\n/g' $miyagiOut|  tr -d [] > $converted
