function brukerRawToImages(folderPath, experimentsFolders, gMaxBruker, filename)
% brukerRawToImages - Transform raw Bruker Images to niftis.
% It concatenates them in the order expressed with the experimentsFolders variable. 
% Calculates the xps structure needed for the md-dmro toolbox.
% INPUTS:
%   folderPath: Folder path of the bruker folder with the images
%   experimentsFolders: vector of the desired folder to get the xps and
%           concatenation done
%   gMaxBruker: Max gradeint strength of the bruker system. Needed for the
%               correct B-tensor calculation. Note: I am looking for a way 
%               to extract  it fromparameters
%   filename: file name for the nifti and xps soutput

% Author: Ricardo Rios
% email:  ricardo.rios@cimat.mx

% ToDo an Optional list to erase volumes

xps = xpsFromBrukerRaw(folderPath, experimentsFolders, gMaxBruker);
allData = extractDataFromBrukerRaw(folderPath, experimentsFolders);

save([filename '_xps'], 'xps')
niftiwrite(allData,filename)

end