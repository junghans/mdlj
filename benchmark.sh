#!/bin/bash

make
cut="1.3"
t0=$(mktemp mdljout.XXXX)
for i in " " "-rvl $cut" "-rlc $cut" "-rvl $cut -rlc $cut"; do
  t1=$(mktemp mdljout.XXXX)
  echo ./mdlj $i "$@"
  time ./mdlj $i "$@" | sed '/^#/d' > "$t1"
  echo "-------------------------"
  [[ -s $t0 ]] || cp "$t1" "$t0"
  cmp "$t0" "$t1" || echo "Results differ!!!"
  rm "$t1"
done
rm "$t0"
