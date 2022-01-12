#!/bin/bash

mode=$1
SAMPLES_MC=( \
    QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8  \
    QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8  \
    QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8  \
    QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8  \
    QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
    QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8  \
    QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8  \
    QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
)

for sample in ${SAMPLES_MC[*]}; do
	./build_objects ntuples/dasntuples_${sample}.root output_bkg_${sample}_${mode}.root 0 0 $mode
done
rm output_allbkg_${mode}.root
hadd output_allbkg_${mode}.root output_bkg_*${mode}.root
rm output_bkg_*${mode}.root
