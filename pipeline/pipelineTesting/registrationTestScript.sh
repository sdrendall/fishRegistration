#!/bin/bash

inputPath='testInput.mhd'
sliceIndex='191.0'
refOutPath='testRefOut.tiff'
annoOutPath='testAnnoOut.mhd'
hemiOutPath='testHemiOut.tiff'
logPath='testLogOut.txt'
atlasBase='/home/sam/Dropbox/grayLab/allenReferenceAtlas_mouseCoronal/atlasVolume'
refAtlasPath="$atlasBase/atlasVolume.mhd"
refAnnoPath="$atlasBase/annotation.mhd"
refHemiPath="$atlasBase/hemisphereMask.mhd"

../registerSliceToAtlas "$inputPath" "$sliceIndex" "$refOutPath" "$annoOutPath" "$hemiOutPath" "$logPath" "$refAtlasPath" "$refAnnoPath" "$refHemiPath"
