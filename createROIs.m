function [] = createROIs(ip)
%% Create event ROIs from raw data
nExperiments = length(ip.allExperiments);
logCellEventTimes = struct;

fprintf('Creating event ROIs...\n');
cd (ip.dataDir);
listindex = 0;
for ee = 1:nExperiments
    cd (ip.allExperiments{ee});
    if ip.roiAll_TS == 1
        dirData = dir; % Get the data for the current directory
        dirData(strncmp({dirData.name}, '.', 1)) = []; % Remove false files
        dirData(strncmp({dirData.name}, '~', 1)) = [];
        dirIndex = [dirData.isdir];  % Find the index for directories
        tsDirs = {dirData(dirIndex).name}';
        tsDirsField = tsDirs;
    elseif ip.roiAll_TS == 0
        dirData = dir; % Get the data for the current directory
        idx = 1;
        for ii = 1:length(dirData) % Include only selected files
            if (dirData(ii).name(length(dirData(ii).name)) == '#') == 1
                tsDirs(idx,1) = {dirData(ii).name};
                tsDirsField(idx,1) = {dirData(ii).name(1:(length(dirData(ii).name)-1))};
                idx = idx+1;
            end
        end
    end
    
    for ts = 1:length(tsDirs)
        listindex = listindex + 1;
        ip.(char(ip.allExperiments{ee})).(char(tsDirsField{ts}))=[];
        listROIs(listindex,1) = ip.allExperiments(ee); %#ok<*AGROW>
        listROIs(listindex,2) = tsDirs(ts);
        listROIs(listindex,3) = tsDirsField(ts);
    end
    cd ../
    clear tsDirs tsDirsField
end

% Import all cellMetadata
listindex = 0;

for ee = 1:nExperiments
    tsDirs = fieldnames(ip.(char(ip.allExperiments{ee})));
    for ts = 1:length(tsDirs)
        listindex = listindex + 1;
        pw=[ip.dataDir filesep ip.allExperiments{ee} filesep listROIs{listindex,2} filesep 'matlabOutput'];
        cd (char(pw));
        position = listROIs{listindex,2};
        if (position(length(position)) == '#') == 1
            position = position(1:(length(position)-1));
        end
        
        filename = char(sprintf('%s%s', position, '_cellDataStructure'));
        load(filename);
        
        experiment = ip.allExperiments{ee};
        nCells = length(fieldnames(cellDataStructure.(experiment).(position)));
        eventOrder = {'first' 'second' 'third' 'fourth' 'fifth' 'sixth' 'seventh' 'eighth' 'nineth' 'tenth'};
        
        for iCells = 1:nCells
            coordinateMatrix = [];
            nEvents = 0;
            if iCells < 10
                cellID = sprintf('cell_00%d',iCells);
            elseif iCells < 100
                cellID = sprintf('cell_0%d',iCells);
            else
                cellID = sprintf('cell_%d',iCells);
            end
            
            nCellEvents = length(fieldnames(cellDataStructure.(experiment).(position).(cellID).allEvents));
            

            for iCellEvents = 1:nCellEvents
                
                eventStatus = cellDataStructure.(experiment).(position).(cellID).allEvents.(eventOrder{iCellEvents}).status;
                eventFrame = cellDataStructure.(experiment).(position).(cellID).allEvents.(eventOrder{iCellEvents}).frame;
                validationStatus = cellDataStructure.(experiment).(position).(cellID).allEvents.(eventOrder{iCellEvents}).validationStatus;
                cellCoordinates = cellDataStructure.(experiment).(position).(cellID).trackingCoordinates;
                
                switch eventStatus
                    case {'ab inito', 'detected'}
                        
                                coordinateMatrix(nEvents+1,1) = floor(cellCoordinates{eventFrame,1});   % x-coordinate
                                coordinateMatrix(nEvents+1,2) = floor(cellCoordinates{eventFrame,2});   % y-coordinate
                                coordinateMatrix(nEvents+1,3) = eventFrame;                             % t-coordinate
                                coordinateMatrix(nEvents+1,4) = iCells;                                 % cell ID
                                coordinateMatrix(nEvents+1,5) = iCellEvents;                            % numbers of previously excluded events in cell
                                nEvents = nEvents + 1;
                end % switch eventstatus
            
            end % iCellEvents
            if isempty(coordinateMatrix) == 0
            ip.(ip.allExperiments{ee}).(tsDirs{ts}).cellMetadata.(cellID).data = coordinateMatrix;
            end
        end % iCells
        
        cd ../../..
    end % tsDirs
    clear tsDirs tsDirsField
end % nExperiments


% Subpartition all cellMetadata to eventData
% Create eventMetadata subdirectories
% Write eventMetadata.mat files to subdirectories
listindex = 0;
showList = cell(0);
for ee = 1:nExperiments
    tsDirs = fieldnames(ip.(char(ip.allExperiments{ee})));
    
    for ts = 1:length(tsDirs)
        eventCounter = 0;
        listindex = listindex + 1;
        showList(length(showList)+1,1) = {sprintf('%s\n', char(listROIs{listindex,2}))};
        % Import background mask data and generate background mask

        cellDirs = fieldnames(ip.(char(ip.allExperiments{ee})).(char(tsDirs{ts})).cellMetadata);
        pw=[ip.dataDir filesep char(ip.allExperiments{ee}) filesep listROIs{listindex,2}];
        cd (char(pw));
        
         destination = 'eventMetadata';
        if ~exist(destination,'dir')
            mkdir(destination);
        end
        cd 'eventMetadata'
        
        nCellEvents = 1;
        allCellFields = fieldnames(ip.(ip.allExperiments{ee}).(tsDirs{ts}).cellMetadata);
        for iCellDirs = 1:length(allCellFields)
            
            cellName = allCellFields{iCellDirs};
            
            mkdir (char(cellName))
            cd (char(cellName))
            
            coordinateMatrix = ip.(char(ip.allExperiments{ee})).(char(tsDirs{ts})).cellMetadata.(cellName).data;
            for tt = 1:size(coordinateMatrix,1)
                
                eventMetadata = coordinateMatrix(tt,:);

                cellEventOrder = coordinateMatrix(tt,5);
                if cellEventOrder < 10
                    eventName = sprintf('%s_event_0%d',cellName,cellEventOrder);
                else
                    eventName = sprintf('%s_event_%d',cellName,cellEventOrder);
                end

                mkdir (char(eventName))
                cd (char(eventName))
                showList(length(showList)+1,1)= {sprintf('-%s',char(eventName))};
                
                save('eventMetadata.mat', 'eventMetadata');
             
                ip.(ip.allExperiments{ee}).(tsDirs{ts}).eventMetadata.(cellName).(eventName).eventMetadata = eventMetadata;
                ip.(ip.allExperiments{ee}).(tsDirs{ts}).eventMetadata.(cellName).(eventName).eventDirPath = pwd;
                logCellEventTimes.(ip.allExperiments{ee}).(tsDirs{ts}).(cellName)(nCellEvents,1) = eventMetadata(1,3);
                logCellEvents.(ip.allExperiments{ee}).(tsDirs{ts}).(cellName).(eventName) = 0;
                ip.(ip.allExperiments{ee}).(tsDirs{ts}).eventCounter = eventCounter + 1;
                
                nCellEvents = nCellEvents + 1;
                eventCounter = eventCounter + 1;
                cd ../
            end
            cd ../
        end
        clear cellDirs
        cd ../../..
    end
    clear tsDirs
end

fprintf('Adding events to queue...\n')
disp(showList)
fprintf('\n')


% Import raw data images as specified by eventMetadata, create and save ROIs
for ee = 1:nExperiments
    fprintf('\nCreating ROIs...\n')
    ip.iExperiments = ee;
    allChannels = ip.allChannels.(char(ip.allExperiments(ee)));
    ip.tsChannels = allChannels;
    nChannels = length(allChannels);
    ip.nChannels = nChannels;
    tsDirs = fieldnames(ip.(char(ip.allExperiments{ee})));
    
    for ts = 1:length(tsDirs)
        if ip.roiAll_TS == 0
            tsDirPath = sprintf('%s#',tsDirs{ts});
        else
            tsDirPath = tsDirs{ts};
        end
        
        % Determine if files should be transferred to SSD drive for processing
        ip.decon_path = [ip.dataDir filesep ip.allExperiments{ee} filesep tsDirPath filesep 'decon'];
        cd (ip.fdp)
        benefit = 'no';
        switch ip.fileAccessMode
            case 'Optimize'
                switch benefit
                    case 'yes'
                        cd(ip.fdp)
                        copymoveSSD('decon', 'deconSSD', ip);
                        pDeconPath = [ip.SSD filesep 'Processing directory' filesep 'resampledImages'];
                    case 'no'
                        pDeconPath = [ip.dataDir filesep ip.allExperiments{ee} filesep tsDirPath filesep 'resampledImages'];
                end
            case 'NAS'
                pDeconPath = [ip.dataDir filesep ip.allExperiments{ee} filesep tsDirPath filesep 'resampledImages'];
        end
        
        
        nCellEvents = 0; %#ok<NASGU>
        cellDirs = fieldnames(ip.(char(ip.allExperiments{ee})).(char(tsDirs{ts})).eventMetadata);
        for iCellDirs = 1:length(cellDirs)
            cellEvents = fieldnames(ip.(char(ip.allExperiments{ee})).(char(tsDirs{ts})).eventMetadata.(char(cellDirs{iCellDirs})));
            
            for iCellEvents = 1:length(cellEvents)
                fprintf('  Processing %s\n', char(cellEvents(iCellEvents)))
                
                eventMetadata = ip.(char(ip.allExperiments{ee})).(char(tsDirs{ts})).eventMetadata.(char(cellDirs{iCellDirs})).(char(cellEvents{iCellEvents})).eventMetadata;
                
                cd (char(pDeconPath));
                
                % Constrain ROI boudaries to within dataset dimensions
                xCC = eventMetadata(:,1);
                yCC = eventMetadata(:,2);
                tCC = eventMetadata(:,3); %%%%%%%%%%%%%%%%
                
                xLim = xCC - ip.extROI;
                xEnd = xCC + ip.extROI;
                yLim = yCC - ip.extROI;
                yEnd = yCC + ip.extROI;
                
                if xLim < 1
                    xLim = 1;
                end
                if xEnd > ip.xDim.(char(ip.allExperiments{ee}))
                    xEnd = ip.xDim.(char(ip.allExperiments{ee}));
                end
                
                if yLim < 1
                    yLim = 1;
                end
                if yEnd > ip.xDim.(char(ip.allExperiments{ee}))
                    yEnd = ip.xDim.(char(ip.allExperiments{ee}));
                end
                
                tStart = tCC - ip.extTimePre;
                tEnd = tCC + ip.extTimePost;
                
                % Start reading and cropping images
                for iChannels = 1:nChannels
                    cd (char(allChannels{iChannels}))
                    
                    clear dir
                    dirData = dir; % Get the data for the current directory
                    dirData(strncmp({dirData.name}, '.', 1)) = []; % Remove false files
                    allFiles = {dirData.name}';
                    
                    % Find number of times from file names
                    scan = 1;
                    iTime = 0;
                    while scan == 1
                        iTime = iTime + 1;
                        if iTime < 10
                            tStr = 't00%d';
                        elseif iTime < 100
                            tStr = 't0%d';
                        else
                            tStr = 't%d';
                        end
                        
                        time = sprintf(tStr, iTime);
                        for iFile = 1:size(allFiles,1)
                            if  isempty(strfind(allFiles{iFile,1},time)) == 0
                                break
                            elseif iFile == size(allFiles,1)
                                scan = 0;
                                nTimes = iTime - 1;
                            end
                        end
                    end
                    
                    % Constrain time boudaries to within dataset dimensions
                    if iChannels == 1
                        if tStart < 1
                            tStart = 1;
                        end
                        if tEnd > nTimes
                            tEnd = nTimes;
                        end
                    end
                    
                    % Find number of planes from file info
                    clear dir
                    dirData = dir; % Get the data for the current directory
                    dirData(strncmp({dirData.name}, '.', 1)) = []; % Remove false files
                    dirData(strncmp({dirData.name}, '~', 1)) = [];
                    allFiles = {dirData.name}';
                    for iFile = 1:size(allFiles,1)
                        if  isempty(strfind(allFiles{iFile},'t001')) == 0
                            break
                        end
                    end
                    nPlanes = size(imfinfo(char(allFiles(iFile))),1);
                    
                    for iTime = tStart:tEnd
                        syntax = '%s';
                        if iTime < 10
                            stringT = 't00%d';
                        elseif iTime < 100
                            stringT = 't0%d';
                        else
                            stringT = 't%d';
                        end
                        timeForm = sprintf(char(syntax), stringT);
                        timeStamp = sprintf(char(timeForm),iTime);
                        for findImage = 1:size(allFiles,1)
                            if isempty(cell2mat(strfind(allFiles(findImage), timeStamp))) == 0
                                image = zeros(ip.xDim.(char(ip.allExperiments{ee})), ip.yDim.(char(ip.allExperiments{ee})), nPlanes);
                                for iPlane = 1:nPlanes
                                    image(:,:,iPlane) = imread(char(allFiles(findImage)), iPlane);
                                end
                                break
                            end
                        end
                        if iTime == tStart && iChannels == 1
                            eventROI = zeros(yEnd-yLim+1, xEnd-xLim+1, nPlanes, (tEnd-tStart), nChannels);
                        end
                        zROI = image(yLim:yEnd, xLim:xEnd,:);
                        eventROI(:,:,:,(iTime-tStart+1), iChannels) = zROI;
                    end
                    cd ../
                end
                
                % create new metadata file for tracking with coordinates
                % adjusted to the ROI crop -> save in eventROI folder also
                % add features with controlled search distances
                
                % Image X-dim = matrix dim(2)
                % -------
                % -------
                % -------
                
                % Image Y-dim = matrix dim(1)
                % |||||||
                % |||||||
                % |||||||
                
                % Save new coordinates of the event (as in ROI crop)
                % Correct xLim, yLim, tStart by subtracting 1
                % M(0,0) does not exist
                roiMetadata = struct;
                roiMetadata.xCC = eventMetadata(1,1) - (xLim-1);    % xCC_roi
                roiMetadata.yCC = eventMetadata(1,2) - (yLim-1);     % yCC_roi
                roiMetadata.tCC = eventMetadata(1,3) - (tStart-1);    % tCC_roi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                roiMetadata.xLim = xLim;
                roiMetadata.xEnd = xEnd;
                roiMetadata.yLim = yLim;
                roiMetadata.yEnd = yEnd;
                roiMetadata.tStart = tStart;
                roiMetadata.tEnd = tEnd;
                roiMetadata.nPlanes = nPlanes;
                roiMetadata.xCC_ORG = xCC;
                roiMetadata.yCC_ORG = yCC;
                roiMetadata.tCC_ORG = tCC; %#ok<*STRNU>
                roiMetadata.allEventFramesIndex = eventMetadata(1,4);
                roiMetadata.prevExcludedEvents = eventMetadata(1,5);
                
                cd (char(pDeconPath))
                cd ../
                
                switch benefit
                    case 'yes'
                        destination = 'eventROIsSSD';
                    case 'no'
                        destination = 'eventROIs';
                end
                
                if ~exist(destination,'dir')
                    mkdir(destination);
                end
                cd    (char(destination))
                
                destination = (char(cellDirs{iCellDirs}));
                if ~exist(destination,'dir')
                    mkdir(destination);
                end
                cd    (char(destination))
                
                mkdir (char(cellEvents{iCellEvents}))
                cd    (char(cellEvents{iCellEvents}))
                
               files = {'eventROI' 'roiMetadata'};
                
                for iFiles = 1:length(files)
                    if ~exist(files{iFiles},'file') == 0
                        delete(files{iFiles})
                    end
                    save(files{iFiles}, files{iFiles})
                end
            end % iCellEvents
        end
        switch benefit
            case 'yes'
                cd (ip.fdp)
                clearSSD('deconSSD', ip)
                ip.eventROIs_path = [ip.dataDir filesep ip.allExperiments{ee} filesep tsDirPath filesep 'eventROIs'];
                copymoveSSD('eventROIsSSD', 'eventROIs', ip);
                clearSSD('eventROIsSSD', ip);
            case 'no'
        end
        
    end
end

cd(ip.fdp)
fprintf('\nROIs have been created.\n');
end


