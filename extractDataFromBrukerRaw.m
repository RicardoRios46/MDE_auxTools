function allData = extractDataFromBrukerRaw(folderPath, exps)
% extractDataFromBrukerRaw - Extract bruwerRaw file into matrixes,
% concatenate them all in the same matrix for easy usage

% Author: Ricardo Rios
% email:  ricardo.rios@cimat.mx

proc_num = 1; % variable to play if another recon is needed

% Start empy dataset
allData = [];

% Loop between experiments
for i_exp = 1:length(exps)  
    expno = exps(i_exp);
    % read Bruker
    source = strcat(folderPath, filesep, num2str(expno),filesep, ...
        'pdata', filesep, num2str(proc_num), filesep);
    visu_params = readBrukerParamFile(strcat(source, 'visu_pars'));
    [imagee, ~] = readBruker2dseq(strcat(source,'2dseq'), visu_params);

    % Orientation (just for matlab plotting)
    permuteOrder=1:ndims(imagee);
    permuteOrder([1 2])=permuteOrder([2 1]);
    shape_imagee = size(imagee);
    data = flip(permute(reshape(imagee,[shape_imagee(1) shape_imagee(2) shape_imagee(3) shape_imagee(5)]),permuteOrder),2);

    allData = cat(4, allData, data);       
end

end