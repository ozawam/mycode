#!/bin/bash
MESSAGE=${1:-"add a subtle touch"}

#git add Makefile git_automatic.sh problem.c
git add *
git commit -m "${MESSAGE}"
#git push origin master
git push origin no_collision_I1992

