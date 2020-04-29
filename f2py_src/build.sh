#!/bin/bash


top=$("pwd")
echo $top
for d in */ ; do
  echo $d
  cd $d
  if test -f "build.sh"; then
    ./build.sh
    cp *.so ../../f2py_lib/
  fi
  cd $top
done
