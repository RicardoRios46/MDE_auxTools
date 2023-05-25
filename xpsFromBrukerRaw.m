function xps = xpsFromBrukerRaw(folderPath, exps, gMaxBruker)
% xpsFromBrukerRaw - Builds the xps structure needed in the MD-dMRI toolbox
% from the bruker folder with the indicated experiments (exps). 
% It needs gMaxBruker for the correct calculation of the B-tensor, I am looking
% for a way to extract it from methods parameters

% Author: Ricardo Rios
% email:  ricardo.rios@cimat.mx

% Intialitating xps structure
xps.n = 0;
xps.b = [];
xps.bt = [];
xps.b_delta=[];
xps.b_eta=[];
xps.u=[];
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
    bTensorVec = zeros(subTotalImages,6);
    
    bValue = zeros(subTotalImages,1); % not in SI units. use b
    b = zeros(subTotalImages,1);
    b_delta = zeros(subTotalImages,1);
    b_eta = zeros(subTotalImages,1);
    u = zeros(subTotalImages,3);

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

        % Calculate B-tensor
        [bTensor, bValue(j_dir)] = calculateBTensor(waveformsConcatenated, rasterTime);
        bTensorVec(j_dir,:) = tm_3x3_to_1x6(bTensor);
        

        % Extract b_tensor parameters
        btensor_params = tm_3x3_to_tpars(bTensor);

        b(j_dir)        = btensor_params.trace;
        b_delta(j_dir)  = btensor_params.delta;
        b_eta(j_dir)    = btensor_params.eta;
        u(j_dir,:) = get_u_from_tenParam(btensor_params);
        


    end
    xps.n = xps.n + subTotalImages;
    xps.b = [xps.b; b];
    xps.bt = [xps.bt; bTensorVec];
    xps.b_delta=[xps.b_delta; b_delta];
    xps.b_eta=[xps.b_eta; b_eta];
    xps.u=[xps.u; u];
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
    % extract acquisition matrix (fov rotations)
    R_acqp = squeeze(method_params.Dw_acqgradmatrix_override(index,:,:));
    % extract rotation matrix for diffusion directions
    R_diff = zeros(3,3);
    R_diff(1,1)=method_params.DwR00(index);
    R_diff(1,2)=method_params.DwR01(index);
    R_diff(1,3)=method_params.DwR02(index);
    R_diff(2,1)=method_params.DwR10(index);
    R_diff(2,2)=method_params.DwR11(index);
    R_diff(2,3)=method_params.DwR12(index);
    R_diff(3,1)=method_params.DwR20(index);
    R_diff(3,2)=method_params.DwR21(index);
    R_diff(3,3)=method_params.DwR22(index);

    R = R_acqp * R_diff; %todo verify this in budde sequence
end

function wfInterpolated = fastInterpolation(wf, duration, rasterTime, dt)
    xq = 0:rasterTime:duration ; % Time points to interpolate
    x =  0:dt:duration ; % time vector 
    wfInterpolated = interp1(x,wf,xq);
end



function u = get_u_from_tenParam(tp)
% code snippet taken from mdm_xps_from_bt function in md-MRI
% toolbox https://github.com/markus-nilsson/md-dmri
    if (tp.eta < 0.1) % skewed tensors?
        
        if (tp.delta == 0) % spherical
            u = tp.lambda33vec;
        elseif (tp.delta > 0) % prolate-to-stick
            u = tp.lambda33vec;
        elseif (tp.delta < 0) % oblate
            u = tp.lambda11vec;
        end
    else
        % not defined for assymmetric tensors
        u = [NaN NaN NaN];
    end
end


function t = tm_3x3_to_1x6(t)
% function t = tm_3x3_to_1x6(t)
% function taken from md-MRI toolbox https://github.com/markus-nilsson/md-dmri
%
% Convert a second-order 3x3 tensor to  Voigt-format 1x6

t = t([1 5 9 2 3 6]) .* [1 1 1 sqrt(2) sqrt(2) sqrt(2)];

end



function t = tm_3x3_to_tpars(t3x3)
% function t = tm_3x3_to_tpars(t3x3)
% function taken from md-MRI toolbox https://github.com/markus-nilsson/md-dmri
%
% Calculate derived tensor parameters.
%
% See definitions in
% Topgaard, In Diffusion NMR of Confined Systems: Fluid Transport in Porous Solids and Heterogeneous Materials;
% Valiullin, R., Ed.; Royal Society of Chemistry: Cambridge, UK, 2016, p 226.
% http://dx.doi.org/10.1039/9781782623779-00226


if (any(isnan(t3x3(:))) || any(isinf(t3x3(:))))
    t3x3 = zeros(3,3);
end

t.t1x6 = tm_3x3_to_1x6(t3x3);

[V,D] = eig(t3x3);
lambdas = diag(D);

% Trace and isotropic average
t.trace = sum(lambdas);
t.iso   = t.trace/3;

[~,ind] = sort(lambdas,'descend');

% Eigenvalue ordering convention as in Eq. (7.2) Topgaard 2016
t.lambda33    = lambdas(ind(1));                        
t.lambda22    = lambdas(ind(2));                       
t.lambda11    = lambdas(ind(3));  
t.lambda33vec = V(:,ind(1))';
t.lambda22vec = V(:,ind(2))';
t.lambda11vec = V(:,ind(3))';

% Linear, planar, and spherical components, see Eq. (7.5) Topgaard 2016
t.l = t.lambda33 - t.lambda22;
t.p = 2*(t.lambda22 - t.lambda11);
t.s = 3*t.lambda11;

% Eigenvalue ordering convention as in Eq. (7.10) Topgaard 2016
Dlambdas = abs(lambdas-t.iso);
[~,ind] = sort(Dlambdas,'descend');
t.lambdazz = lambdas(ind(1));
t.lambdaxx = lambdas(ind(2));                        
t.lambdayy = lambdas(ind(3));
t.lambdazzvec = V(:,ind(1))';
t.lambdaxxvec = V(:,ind(2))';
t.lambdayyvec = V(:,ind(3))';

% Normalized anisotropy and asymmetry , see Eq. (7.11) Topgaard 2016
t.delta = (t.lambdazz - (t.lambdayy+t.lambdaxx)/2)./(t.trace + eps);
t.eta = 3*(t.lambdayy - t.lambdaxx ) ./ (2*t.trace.*t.delta+eps);

% Variance of eigenvalues
t.vlambda = 1/3*((t.lambda33 - t.iso).^2 + (t.lambda22 - t.iso).^2 + (t.lambda11 - t.iso).^2);

% Fractional anisotropy, see Basser. J. Magn. Reson. B 111, 209 (1996).
t.fa = sqrt(1/2)*sqrt((t.lambda33-t.lambda22).^2+(t.lambda33-t.lambda11).^2+(t.lambda22-t.lambda11).^2)...
                    ./sqrt(t.lambda33.^2+t.lambda22.^2+t.lambda11.^2);
t.fa(isnan(t.fa)) = 0;

% Westin's macroscopic anisotropy measure, see Westin. Neuroimage 135, 345 (2016).
t.cm = t.vlambda./t.iso.^2;
t.cm(isnan(t.cm)) = 0;

% Westin's linear, planar, spherical measures, see Westin. Med. Image Anal. 6, 93 (2002).
t.cl = t.l./t.trace;
t.cl(isnan(t.cl)) = 0;
t.cp = t.p./t.trace;
t.cp(isnan(t.cp)) = 0;
t.cs = t.s./t.trace;
t.cs(isnan(t.cs)) = 0;


% Compute parallel and perpendicular components
% First determine whether it is prolate or oblate
if (t.cl > t.cp)
    t.para = t.lambda33;
    t.perp = (t.lambda11 + t.lambda22) / 2;
    t.norm = t.lambda33vec;
else
    t.para = t.lambda11;
    t.perp = (t.lambda22 + t.lambda33) / 2;
    t.norm = t.lambda11vec;
end

end