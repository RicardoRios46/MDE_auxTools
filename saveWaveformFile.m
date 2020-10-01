function saveWaveformFile(Waveforms, gMaxBruker, directionVector, fileName) 
% saveWaveformFile - Save waveform structure in dormatted files for the MCW
% pulse sequences toolbox
waveforms{1} = Waveforms.waveform1'; % Takings transpose for easy fprint
waveforms{2} = Waveforms.waveform2';
duration(1) = Waveforms.waveform1_duration;
duration(2) =  Waveforms.waveform2_duration;
nPoints(1) = Waveforms.waveform1_Npoints;
nPoints(2) = Waveforms.waveform2_Npoints;

for iWf = 1:2
    
    waveformPercentage = waveforms{iWf}./gMaxBruker; % transform into a percentage
    
    fileID = fopen([fileName num2str(iWf)],'w');

    % contant lines needed in the bruker format
    fprintf(fileID,'##TITLE= Waveform\n');
    fprintf(fileID,'##JCAMP-DX= 5.00 Bruker JCAMP library\n');
    fprintf(fileID,'##DATA TYPE= Shape Data\n');
    fprintf(fileID,'##ORIGIN= Bruker Analytik GmbH\n');
    fprintf(fileID,'##OWNER= <nmrsu>\n');
    fprintf(fileID,'##DATE= xx\n');
    fprintf(fileID,'##TIME= xx\n');
    fprintf(fileID,'##MINX= 0\n');
    fprintf(fileID,'##MAXX= 1\n');
    fprintf(fileID,'##MINY= 0\n');
    fprintf(fileID,'##MAXY= 1\n');
    fprintf(fileID,'##$SHAPE_EXMODE= Gradient\n');
    fprintf(fileID,'##$SHAPE_TOTROT= 0\n');
    fprintf(fileID,'##$SHAPE_BWFAC= 0\n');
    fprintf(fileID,'##$SHAPE_INTEGFAC= 0\n');
    fprintf(fileID,'##$SHAPE_MODE= 0\n');

    % print file
    fprintf(fileID,'##DIRECTIONVEC= %d %d %d\n',directionVector);

    fprintf(fileID,'##DURATION= %f\n',duration(iWf));
    fprintf(fileID,'##NPOINTS= %d\n',nPoints(iWf));

    fprintf(fileID,'##XYDATA= (T X Y Z)\n');

    fprintf(fileID,'%f %f %f\n',waveformPercentage);

    fprintf(fileID,'##END');

    fclose(fileID);
end
end