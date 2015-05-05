#!/bin/bash

inputPath='testInput.mhd'
sliceIndex='365.0'
refOutPath='testRefOut.mhd'
annoOutPath='testAnnoOut.mhd'
hemiOutPath='testHemiOut.mhd'
logPath='testLogOut.txt'
atlasBase='/home/sam/code/fishRegistration/atlasVolume'
refAtlasPath="$atlasBase/atlasVolume.mhd"
refAnnoPath="$atlasBase/annotation.mhd"
refHemiPath="$atlasBase/hemisphereMask.mhd"

../registerSliceToAtlas "$inputPath" "$sliceIndex" "$refOutPath" "$annoOutPath" "$hemiOutPath" "$logPath" "$refAtlasPath" "$refAnnoPath" "$refHemiPath"
