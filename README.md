fishRegistration
================

Scripts and MATLAB apps to facilitate the use of the ITK-based image registration software found [here](https://github.com/sdrendall/image_registration_tools)

Setup
-----

These tools require MATLAB and Python 2.7. I used [anaconda python](https://www.anaconda.com/download).

### python dependencies

Multithreading is provided using [`twisted`](https://github.com/twisted/twisted). An optional .xls/.xlsx => tsv/csv converter is included that is based on `xlrd`. Everything else should be included with anaconda.

You also must add the `pythonMods` directory to your `PYTHONPATH`.

```
pip install twisted xlrd
```

### matlab

These scripts make use of several third party matlab toolboxes which are included in `matlabToolboxes` for convenience. For these scripts to function, you must add the `matlabFcns` and `matlabToolboxes` directories to your matlab path.

### image\_registration\_tools

You also must compile and install [ITK](https://github.com/InsightSoftwareConsortium/ITK) and [the image\_registration\_tools](https://github.com/sdrendall/image_registration_tools). The image\_registration\_tools binaries must be somewhere on your `PATH` for the registerImages script to work.

### brain atlases

We use a modified version of the [allen brain atlas](http://www.brain-map.org) We also use a hemisphere mask for lateralization of allen brain annotations (which are symmectrical in the official allen reference atlas). You can download them [here](https://www.dropbox.com/s/lwxk3p0r5mbju8q/atlasVolume.zip?dl=0)

Usage
-----

These tools create a `.registrationData` directory where intermediate and output images/data are stored throughout the registration process. There are three major steps to the process:

1. Preprocessing -- `preprocessImages` -- Images are converted from their original file format to `.mhd`, the pixel intensities are normalized, the background (the area outside of the brain section) is zeroed, and they are scaled down to the size of the atlas (320 x 456)
2. Coordinate Indexing -- `assignBregmaCoordinates` or `importAtlasIndicies.py` -- Images are labelled with position along the rostral/caudal axis.
3. Registration -- `registerImages` -- Images are registered to the corresponding of the allen brain atlas.
