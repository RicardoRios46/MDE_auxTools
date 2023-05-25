# MDE_auxTools (Beta)

*MDE_auxTools* contains some useful functions that help to implement B-tensor Encoding sequences on a Bruker scanner. The functions here are intended to use as connectors for other repositories developed by the DW-MRI scientific community to have a complete acquisition pipeline:

![image](https://drive.google.com/uc?export=view&id=1xVbbHms5eVItVbdXMG5nyviWYzp6GRrk)

Use *help* in Matlab command window to learn how to use the functions and the inputs needed.
## Repositories used in the pipeline:
* [NOW](https://github.com/jsjol/NOW)
* [MCW sequences](http://osf.io/ngu4a/)
* [BrkRaw](https://github.com/BrkRaw/bruker)
* [MD-dMRI](https://github.com/markus-nilsson/md-dmri)

## Requirements
*nowToSequence* needs that *NOW* toolbox is already on matlab path.



## Notes
*brukerRawToImages* is still a prototype. Ideally takes raw bruker data and outputs a nifti file with its xps for processing. It uses the *pvmatlab* toolbox from Bruker (To get this toolbox you need to contact Bruker support). For now it only works with one shape and one gradient amplitude for each acquisition (it works with multiple gradient rotations).

### Roadmap for future implementations:
* *brukerRawToImages* should work with multiple shapes and amplitudes in each acquisition and use open source software.
* Erase *pvtools* dependency from *extractDataFromBrukerRaw*.
* Optional basic preprocessing steps.
* Tools for basic visualization/evaluation of data and waveforms.
