#!/bin/bash

make
[[ $1 = -v ]] && shift && verbose=yes
cut="1.3"
t0=$(mktemp mdljout.XXXX)
for i in "" "-rvl $cut" "-rlc $cut" "-rvl $cut -rlc $cut"; do
  t1=$(mktemp mdljout.XXXX)
  [[ -z ${i0} ]] && i0="$i"
  echo ./mdlj $i "$@"
  time ./mdlj $i "$@" > "$t1"
  [[ $verbose ]] && cat "${t1}"
  sed -i '/^#/d' "${t1}"
  echo "-------------------------"
  [[ ! -s $t0 ]] && cp "$t1" "$t0" && gold=$i
  if cmp "$t0" "$t1"; then
    [[ $i != $gold ]] && echo "Result of '$i' match with gold (${gold:-no args})"
  else
    echo "Result of '$i' differ from gold (${gold:-no args})" >&2
    exit 1
  fi
  rm "$t1"
done
rm "$t0"
