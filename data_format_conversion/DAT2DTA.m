function DAT2DTA()
% This function converts .DAT file to .DTA and .DTF files

    % Read .DAT file
    [fileName, pathName] = uigetfile('*.DAT','File Selector');
    if ~isequal(fileName,0)
        filePath = strcat(pathName,fileName);
        P = dlmread(filePath);
        x = P(:,1) .* 1e3;
        y = P(:,2);
    
        % Create Brucker files
        TitleString = fileName(1:(end-4));
        filePath = strcat(pathName,TitleString);
        eprsave(filePath,x,y,TitleString)       
    end

end