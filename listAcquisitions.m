function [ processAcquisitions ] = listAcquisitions(ip)
%% Create list of all acquisitions to be processed
nExperiments = numel(ip.allExperiments);

nFiles = 0;
for ee = 1:nExperiments
    pw = [ip.dataDir filesep ip.allExperiments{ee}];
    cd (char(pw))
    if ip.processAll_TS.(char(ip.allExperiments{ee})) == 1
       dirData = dir; % Get the data for the current directory
       dirData(strncmp({dirData.name}, '.', 1)) = []; % Remove false files
       dirData(strncmp({dirData.name}, '~', 1)) = [];
       dirIndex = [dirData.isdir];  % Find the index for directories
       tsDirs = {dirData(dirIndex).name}';   
    elseif ip.processAll_TS.(char(ip.allExperiments{ee})) == 0
       dirData = dir; % Get the data for the current directory
       idx = 1;
       for ii = 1:length(dirData) % Include only selected files
           if (dirData(ii).name(length(dirData(ii).name)) == '#') == 1 
               tsDirs(idx,1) = {dirData(ii).name};
               idx = idx+1;
           end
       end
    end
    
    for ts = 1:length(tsDirs)
        pw = [ip.dataDir filesep ip.allExperiments{ee} filesep tsDirs{ts}];
        cd (char(pw));
        
        clear dir
        dirData = dir; % Get the data for the current directory
        dirData(strncmp({dirData.name}, '.', 1)) = []; % Remove false files
        dirData(strncmp({dirData.name}, '~', 1)) = [];
        
        if exist('offlineStorageTag.mat','file') == 0
            for ii = 1:length(dirData)
                if isempty(strfind(dirData(ii).name, 'Thumbs')) == 1
                    zenFile = dirData(ii).name;
                end
            end
            
            nFiles = nFiles + 1;
            processAcquisitions(nFiles,1) = ip.allExperiments(ee); %#ok<*AGROW>
            processAcquisitions(nFiles,2) = tsDirs(ts);
            processAcquisitions(nFiles,3) = {zenFile};
            processAcquisitions(nFiles,4) = {pwd};
            processAcquisitions(nFiles,6) = {'online'};

            
        else
            nFiles = nFiles + 1;
            processAcquisitions(nFiles,1) = ip.allExperiments(ee); %#ok<*AGROW>
            processAcquisitions(nFiles,2) = tsDirs(ts);
            processAcquisitions(nFiles,3) = {'N/A'};
            processAcquisitions(nFiles,4) = {pwd};
            processAcquisitions(nFiles,6) = {'offline'};
            processAcquisitions(nFiles,7) = {0};
        end
            
    end
    clear tsDirs
end

cd(ip.fdp)
end

