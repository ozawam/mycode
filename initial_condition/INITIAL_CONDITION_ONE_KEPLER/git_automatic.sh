#!/bin/bash
MESSAGE=${1:-"add a subtle touch"}

git add Makefile git_automatic.sh problem.c
git commit -m "${MESSAGE}"
#git push original master
git push original gas_drag
