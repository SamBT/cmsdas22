#!/bin/bash

c++ -lm -o build_objects build_objects.cpp -I include/ `root-config --glibs --cflags`
