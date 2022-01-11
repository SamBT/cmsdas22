#!/bin/bash

mode=$1

./build_objects dasntuples_gg_HH_bbbb_SM.root output_ggHH_SM_${mode}.root 0 1 $mode
