function xps = xpsFromBrukerRaw(folderPath, exps, gMaxBruker)
% xpsFromBrukerRaw - Builds the xps structure needed in the MD-dMRI toolbox
% from the bruker folder with the indicated experiments. 
% It need gMaxBruker for the correct calculation of the B-tensos, I am looking
% for a way to extract it from mehtos parameters

% Author: Ricardo Rios
% email:  ricardo.rios@cimat.mx

% Intialitating xps structure
xps.n = 0;
xps.b = [];
xps.bt = [];
xps.b_delta=[];
xps.s_ind = [];

% Loop between experiments
for i_exp = 1:length(exps)  
    expno = exps(i_exp);
%     % Define paths
%     MethodFilePath = [folder_name, filesep, num2str(expno), ...
%         filesep,'method'];
%     % read Bruker
%     source = strcat(folder_name, filesep, num2str(expno),filesep, ...
%         'pdata', filesep, num2str(proc_num), filesep);
%     acq_params = readBrukerParamFile(strcat(folder_name, filesep,...
%         num2str(expno),filesep,'acqp'));
    method_params = readBrukerParamFile(strcat(folderPath, filesep, ...
        num2str(expno),filesep,'method')); % all the sequence parameters should be here
    
    % Extracting params
    % ToDo if you are looping over shapes dims change
    waveform1= method_params.DwGradShapeArray1; 
    waveform2= method_params.DwGradShapeArray2; 

    wf1Npoints=method_params.DwShapePoints1;
    wf2Npoints=method_params.DwShapePoints2;

    %wf1= squeeze(method_params.DwGradShapeArray1(1,:,:)); % ToDo are we sure about the shape dim?
    %wf2= squeeze(method_params.DwGradShapeArray2(1,:,:)); % ToDo are we sure about the shape dim?

    scaleFactor=method_params.DwGradAmpScale/100;

    % loop over directions
    subTotalImages = method_params.DwNB0 + method_params.DwNDirs;
    bValue = zeros(subTotalImages,1);
    bTensorVec = zeros(subTotalImages,6);
    b_delta = zeros(subTotalImages,1);
    for j_dir = 1:subTotalImages
        % extract Rotation matrix used
        R = getR(method_params, j_dir);

        % rotate waveforms
        waveform1rotated = rotateWaveform(waveform1, R);
        waveform2rotated = rotateWaveform(waveform2, R);

        % Transform from percentage to T/m and scale
        waveform1scaled = waveform1rotated * gMaxBruker * scaleFactor;
        waveform2scaled = waveform2rotated * gMaxBruker * scaleFactor;

        % Transform to s and obtain time_step
        % ToDo check that this correspond to raster time used in other function
        dt1 = method_params.DwGradDur1 * 1e-3 / (wf1Npoints(1) -1);
        dt2 = method_params.DwGradDur2 * 1e-3 / (wf2Npoints (1)-1);

        % ideally dt1 and dt2s should be of the same.
        % interpolating waveforms if thats not the case
        % ToDo check this is actually working    
        if dt1 == dt2
            rasterTime = dt1;
        else
            rasterTime = 5e-6;
            waveform1scaled = fastInterpolation(waveform1scaled, method_params.DwGradDur1 * 1e-3, rasterTime, dt1);
            waveform2scaled = fastInterpolation(waveform2scaled, method_params.DwGradDur2 * 1e-3, rasterTime, dt2);
        end

        % concatenate both waveforms and blank space for RF pulse,
        % also inverting back second waveform
        nZeroSeparationTime = round(method_params.DwGradTsep*1e-3 / rasterTime);
        waveformsConcatenated = [waveform1scaled; zeros(nZeroSeparationTime,3); waveform2scaled*(-1)];

        % Calculate B
        [bTensor, bValue(j_dir)] = calculateBTensor(waveformsConcatenated, rasterTime);
        bTensorVec(j_dir,:) = tensor2vector(bTensor);
        b_delta(j_dir) = calculateBDelta(bTensor);
    end
    xps.n = xps.n + subTotalImages;
    xps.b = [xps.b; bValue];
    xps.bt = [xps.bt; bTensorVec];
    xps.b_delta=[xps.b_delta; b_delta];
    xps.s_ind = [xps.s_ind; i_exp];

end

end

function g = rotateWaveform(waveform, R)
    n = size(waveform,1);
    g = zeros(size(waveform));
    for k = 1:n
        g(k,:) = R*waveform(k,:)';
    end
end

function R = getR(method_params, index)
    R = zeros(3,3);
    R(1,1)=method_params.DwR00(index);
    R(1,2)=method_params.DwR01(index);
    R(1,3)=method_params.DwR02(index);
    R(2,1)=method_params.DwR10(index);
    R(2,2)=method_params.DwR11(index);
    R(2,3)=method_params.DwR12(index);
    R(3,1)=method_params.DwR20(index);
    R(3,2)=method_params.DwR21(index);
    R(3,3)=method_params.DwR22(index);
end

function wfInterpolated = fastInterpolation(wf, duration, rasterTime, dt)
    xq = 0:rasterTime:duration ; % Time points to interpolate
    x =  0:dt:duration ; % time vector 
    wfInterpolated = interp1(x,wf,xq);
end

function vec = tensor2vector(tensor)
    vec = zeros(1,6);
    vec(1) = tensor(1,1);
    vec(2) = tensor(2,2);
    vec(3) = tensor(3,3);
    vec(4) = tensor(1,2);
    vec(5) = tensor(1,3);
    vec(6) = tensor(2,3);
end

function bDelta = calculateBDelta(B)
    bValue = trace(B);
    eigen = eig(B);
    [lambdaAxial, lambdaAxialIndex] = max(eigen);
    eigen(lambdaAxialIndex)=0;
    lamndaRadial = sum(eigen)/2;
    if bValue == 0
        bDelta = 0; % ToDo what is bdelta when bvalue is 0
    else
        bDelta = (lambdaAxial - lamndaRadial)/bValue;
    end
end