function [Waveforms] = nowToSequence(problemNOW, resultNOW, param)
% nowToSequence - Prints usable gradients files for MCW pulse sequences
% toolbox from NOW objects from the NOW toolbox.
% INPUTS:
% problemNOW: problem optimization object. It comes from a NOW optimization
% resultNOW: result structure. It comes from a NOW optimization
% param: Structure with different input parameters for the pipeline:
%   gMaxBruker - Max gradient strength on bruker sistem (T/m)
%   rasterTime - Min raster time on bruker sistem- recommended 10us (s)
%   desiredB - Desired b value for the output waveform (s/mm^2). It must not exceed
%           gradient 100 strength, a warning is displayed if that is the case
%   tol - Tolerance for calculating the desired b value. 10 should do it.
%   directionVector - main direction used in the waveform. Needed if you
%           plan to use the rotations on the MCW sequence toolbox. Check their wiki
%   fileName - Filename of the output files

% Author: Ricardo Rios
% email:  ricardo.rios@cimat.mx

% Scale waveform to desireed b value
scaleFactor =  calculateScaleFactor(resultNOW.g*1e-3, resultNOW.dt, param.desiredB, param.tol);

% Warnings if waveforms use too much gradients
checkGradientUse(resultNOW.g*1e-3,scaleFactor,param.gMaxBruker)

% obtain waveform struct
Waveforms=nowObjectsToWaveforms(problemNOW, resultNOW, param.rasterTime, scaleFactor);

% save to file waveforms
saveWaveformFile(Waveforms, param.gMaxBruker, param.directionVector, param.fileName)

end

