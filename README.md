# MDE_auxTools (Beta)

*MDE_auxTools* contains some useful functions that help to implement B-tensor Encoding sequences on a Bruker scanner. The functions here are intended to use as connectors for other repositories developed by the DW-MRI scientific community to have a complete acquisition pipeline:
![resources](https://docs.google.com/drawings/d/e/2PACX-1vSElx1dUpIEyZj7Qmah8D8nCggIHCiVN6n1rGYA-g6wKyGCstxI22sfRcmqYiMsSHgVOmRbNBAp5AK_/pub?w=734&h=440)

Use *help* in Matlab command window to learn how to use the functions and the inputs needed.
## Repositories used in the pipeline:
* [NOW](https://github.com/jsjol/NOW)
* [MCW sequences](http://osf.io/ngu4a/)
* [MD-dMRI](https://github.com/markus-nilsson/md-dmri)

## Requirements
*nowToSequence* needs that *NOW* toolbox is already on matlab path.

*brukerRawToImages* works with *pvmatlab* toolbox from Bruker. To get this toolbox you need to contact Bruker support.

## Limitations
*MDE_auxTools* is still a work in progress. At this moment *brukerRawToImages* only works if you have used one shape and one gradient amplitude in the MCW sequence toolbox for each acquisition. It works with multiple gradient rotations.

### Roadmap for future implementations:
* *brukerRawToImages* should work with multiple shapes and amplitudes in each acquisition.
* A python version of *brukerRawToImages* that uses [BrkRaw](https://github.com/BrkRaw/bruker) insted of Bruker's *pvmatlab*.
* Function to erase selected volumes from the final Nifti and xps structure.
* Optional basic preprocessing steps.
* Tools for basic visualization/evaluation of data and waveforms.
