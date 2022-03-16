function [] = copyImages(ip)
%% Copy the images to a new folder
experimentFolder = [ip.dataDir filesep ip.listAcquisitions{ip.iAcquisitions,1}];

% Sort raw images into folders
fprintf('\nCopying files to new folder...\n')

readDir = [experimentFolder filesep ip.listAcquisitions{ip.iAcquisitions,2} filesep 'rawTiffs'];
cd (readDir)
clear dir
dirData = dir; % Get the data for the current directory
dirData(strncmp({dirData.name}, '.', 1)) = []; % Remove false files
dirData(strncmp({dirData.name}, '~', 1)) = [];
allFiles = {dirData.name}';

for ii = 1:length(allFiles)
    if contains(allFiles(ii),'c1') == 0
        source = [readDir filesep char(allFiles(ii))];
        destination = [experimentFolder filesep char(ip.listAcquisitions{ip.iAcquisitions,2}) filesep 'readImages' filesep char(allFiles(ii))];
        copyfile(source, destination)
    end
end

end

