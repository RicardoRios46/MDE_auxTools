function Waveforms=nowObjectsToWaveforms(problemNOW, resultNOW, rasterTime, scaleFactor)
%NOW_to_Waveforms - Converts NOW objects into a waveform structure.
%Waveform structue contains the waveforms scaled to scaleFactor, interpolated 
%to rasterTime and separated into two variables for printing in correct 
%format for MCW sequences

% Author: Ricardo Rios
% email:  ricardo.rios@cimat.mx


% ToDo maybe this can be an optional parameter
invertWaveform2 = true;

% scale waveform
waveformScaled = resultNOW.g*scaleFactor; 
waveformScaled = waveformScaled*1e-3; % from mT/m to T/m

% interpolate waveform
waveformScaledInterpolated = interpolateWaveform(waveformScaled, problemNOW, resultNOW, rasterTime);

% extract waveforms
[waveform1, waveform2] = extractWaveforms(waveformScaledInterpolated, problemNOW, rasterTime);
% invert second waveform fo the RF pulse
if invertWaveform2
    waveform2 = waveform2 * (-1); %
end

% output structure
Waveforms.NOW_waveform = resultNOW.g;

Waveforms.waveform1 = waveform1;
Waveforms.waveform2 = waveform2;
Waveforms.waveform1_duration = problemNOW.durationFirstPartActual*1e-3;
Waveforms.waveform2_duration = problemNOW.durationSecondPartActual*1e-3;
Waveforms.waveform1_Npoints = size(waveform1,1);
Waveforms.waveform2_Npoints = size(waveform2,1);

Waveforms.raster_time = rasterTime;

[bTensor, bValue] = calculateBTensor(waveformScaledInterpolated, rasterTime);
Waveforms.bTensor = bTensor;
Waveforms.bValue = bValue;

end

function waveformInterpolated = interpolateWaveform(waveform, problemNOW, resultNOW, rasterTime)
% Simple linear interpolation for the waveform
x =  0:resultNOW.dt:problemNOW.totalTimeActual*1e-3 ; % time vector from NOW problem
xq = 0:rasterTime:problemNOW.totalTimeActual*1e-3 ; % Time points to interpolate

waveformInterpolated = interp1(x,waveform,xq);
end

function [waveform1, waveform2] = extractWaveforms(waveform, problem, rasterTime)
% Extract pair of waveforms (for both sides of the RF pulse) from single waveform
% Waveform in in T/m
%  problem variable can be changed for the actual waveforms durations 

waveform1end = problem.durationFirstPartActual*1e-3 / rasterTime + 1; % Add the last 0
waveform1end = round(waveform1end); % Avoing a real valued index
waveform1 = waveform(1:waveform1end,:);

waveform2start = (problem.totalTimeActual*1e-3 - ...
    problem.durationSecondPartActual*1e-3) / rasterTime + 1; % Extract 1st 0
waveform2start = round(waveform2start); % Avoing a real valued index
waveform2 = waveform(waveform2start:end,:);

end
