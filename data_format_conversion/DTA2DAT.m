function DTA2DAT()
% This function converts .DTA file to .DAT file

    % Read .DTA file
    [fileName, pathName] = uigetfile('*.DTA','File Selector');
    if ~isequal(fileName,0)
        filePath = strcat(pathName,fileName);
        [x, y, Pars] = eprload(filePath)
    
        % Create .DAT file
        fileName = fileName(1:(end-3));
        fileName = strcat(fileName,'dat');
        filePath = strcat(pathName,fileName);
        File = fopen(filePath,'w');
        for i = 1:size(x,1)
            fprintf(File,'%f %f %f\n',x(i),real(y(i)),imag(y(i)));
        end
        fclose(File);
    end

end




 






