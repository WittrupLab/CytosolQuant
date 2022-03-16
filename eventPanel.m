function [panel, panelMetadata] = eventPanel(ip)
%%  Compose time series panel with all events and tracks of the capture
fprintf('Creating event panel...\n');
nChannels = length(ip.allChannels.(char(ip.allExperiments(1))));                  
nExperiments = length(ip.allExperiments);
cumSizeX = ip.spacersy; %#ok<NASGU>
panelfilled = 0;
dMaxPlanes = 0;
dMaxTimes = 0;
dMax_tCC = 0;
dMax_yDim = 0; %#ok<NASGU>
nROIs = 0;
panelMetadata = cell(0);

cd (char(ip.dataDir));

% Create allROI metadata struct
for ee = 1:nExperiments
    pw = [ip.dataDir filesep ip.allExperiments{ee}];
    cd (char(pw));
    if ip.plotAll_TS == 1
       dirData = dir; % Get the data for the current directory
       dirData(strncmp({dirData.name}, '.', 1)) = []; % Remove false files
       dirData(strncmp({dirData.name}, '~', 1)) = [];
       dirIndex = [dirData.isdir];  % Find the index for directories
       tsDirs = {dirData(dirIndex).name}';   
       tsDirsField = tsDirs;
    else
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
        loadCellDataStructure = 1;
        pw = [ip.dataDir filesep ip.allExperiments{ee} filesep tsDirs{ts}];
        cd (char(pw));
        cd 'eventMetadata'
        clear dir
        if ip.plotAll_cells == 1
            dirData = dir; % Get the data for the current directory
            dirData(strncmp({dirData.name}, '.', 1)) = []; % Remove false files
            dirData(strncmp({dirData.name}, '~', 1)) = [];
            dirIndex = [dirData.isdir];  % Find the index for directories
            cellDirs = {dirData(dirIndex).name}';
            cellDirsField = cellDirs;
        else
            dirData = dir; % Get the data for the current directory
            idx = 1;
            for ii = 1:length(dirData) % Include only selected files
                if (dirData(ii).name(length(dirData(ii).name)) == '#') == 1
                    cellDirs(idx,1) = {dirData(ii).name};
                    cellDirsField(idx,1) = {dirData(ii).name(1:(length(dirData(ii).name)-1))};
                    idx = idx+1;
                end
            end
        end
        if exist('cellDirs','var') == 0
                cd ../
                continue
        end
            
        for iCellDirs = 1:length(cellDirs)
            pw = [ip.dataDir filesep ip.allExperiments{ee} filesep tsDirs{ts} filesep 'eventMetadata' filesep cellDirs{iCellDirs}];
            cd (char(pw));
            clear dir
            if ip.plotAll_events == 1
               dirData = dir; % Get the data for the current directory
               dirData(strncmp({dirData.name}, '.', 1)) = []; % Remove false files
               dirData(strncmp({dirData.name}, '~', 1)) = [];
               dirIndex = [dirData.isdir];  % Find the index for directories
               cellEvents = {dirData(dirIndex).name}';
               cellEventsField = cellEvents;
            else
               dirData = dir; % Get the data for the current directory
                idx = 1;
                for ii = 1:length(dirData) % Include only selected files
                    if (dirData(ii).name(length(dirData(ii).name)) == '#') == 1
                        cellEvents(idx,1) = {dirData(ii).name};
                        cellEventsField(idx,1) = {dirData(ii).name(1:(length(dirData(ii).name)-1))};
                        idx = idx+1;
                    end
                end
            end
            if exist('cellEvents','var') == 0
                cd ../
                continue
            end
            
            eventOrder = {'first' 'second' 'third' 'fourth' 'fifth' 'sixth' 'seventh' 'eighth' 'nineth' 'tenth'};
            stopSearch = 0;
            presentEventOrder = 1;
            continueMain = 0;
            for iCellEvents = 1:length(cellEvents)
                
                if stopSearch == 1
                    break
                end
                
                if loadCellDataStructure == 1
                    pw = [ip.dataDir filesep ip.allExperiments{ee} filesep tsDirs{ts} filesep 'matlabOutput'];
                    cd (char(pw));
                    filename = char(sprintf('%s%s', tsDirsField{ts},'_cellDataStructure'));
                    load(filename)
                    loadCellDataStructure = 0;
                end
                
                switch cellDataStructure.(ip.allExperiments{ee}).(tsDirsField{ts}).(cellDirs{iCellDirs}).disqualifyCell
                    case 'yes'
                        break
                end
                
                validationStatus = cellDataStructure.(ip.allExperiments{ee}).(tsDirsField{ts}).(cellDirsField{iCellDirs}).allEvents.(eventOrder{iCellEvents}).validationStatus;

                switch ip.includeEventStatus
                    case 'not validated'
                        switch validationStatus

                            case 'not validated'
                                switch cellDataStructure.(ip.allExperiments{ee}).(tsDirsField{ts}).(cellDirsField{iCellDirs}).eGFPExpressionLevel  
                                    case 'subthreshold'
                                        continue
                                end
                                switch cellDataStructure.(ip.allExperiments{ee}).(tsDirsField{ts}).(cellDirsField{iCellDirs}).controlPosition  
                                    case 'yes'
                                        continue
                                end
                                stopSearch = 1;

                            case 'valid event'
                                if strcmp(eventOrder{presentEventOrder}, ip.includeEventOrder) == 1
                                  stopSearch = 1;%skip event
                                  continue
                                else
                                    presentEventOrder = presentEventOrder+1;
                                    continue
                                end
    
                            case 'false event'
                                
                                if length(cellEvents) < 2 || iCellEvents == length(cellEvents)
                                    stopSearch = 1;%skip event
                                  continue
                                end
                                
                                 for findFirstTrueEvent = 2:length(cellEvents)
                                    scanEvalStatus = cellDataStructure.(ip.allExperiments{ee}).(tsDirsField{ts}).(cellDirsField{iCellDirs}).allEvents.(eventOrder{findFirstTrueEvent}).validationStatus;
                                    switch scanEvalStatus
                                        case 'not validated'
                                            iCellEvents = findFirstTrueEvent;
                                            stopSearch = 1;
                                            break
                                            
                                        case 'false event'
                                            if findFirstTrueEvent == length(cellEvents)
                                                continueMain = 1;
                                                stopSearch = 1;%skip event
                                                break
                                            end
                                            
                                            continue
                                            
                                        case 'valid event'
                                            if strcmp(eventOrder{presentEventOrder}, ip.includeEventOrder) == 1
                                                stopSearch = 1;%skip event
                                                continueMain = 1;
                                                break
                                            elseif findFirstTrueEvent == length(cellEvents)
                                                continueMain = 1;
                                                stopSearch = 1;%skip event
                                                break
                                            else
                                                presentEventOrder = presentEventOrder+1;
                                                continue
                                            end
                                    end
                                 end
                                 
                                if continueMain == 1
                                    continue
                                end
                                
                        end
                        
                    case 'valid events'
                                switch cellDataStructure.(ip.allExperiments{ee}).(tsDirsField{ts}).(cellDirsField{iCellDirs}).eGFPExpressionLevel  
                                    case 'subthreshold'
                                        continue
                                end
                                switch cellDataStructure.(ip.allExperiments{ee}).(tsDirsField{ts}).(cellDirsField{iCellDirs}).controlPosition  
                                    case 'yes'
                                        continue
                                end
                        
                        switch validationStatus
                            case {'not validated', 'false event'}
                                continue
                            case 'valid event'    
                                if strcmp(eventOrder{presentEventOrder}, ip.includeEventOrder) == 1
                                    stopSearch = 1;%add event
                                else
                                    presentEventOrder = presentEventOrder+1;
                                    continue
                                end
                        end
                        
                    case 'false events'
                        switch validationStatus
                            case {'valid event', 'not validated'}
                                continue
                            case 'false event'
                                if presentEventOrder == ip.includeEventOrder
                                    stopSearch = 1;%add event
                                else
                                    presentEventOrder = presentEventOrder+1;
                                    continue
                                end
                        end
                        
                    case 'all events'
                        if presentEventOrder == ip.includeEventOrder
                            stopSearch = 1;%add event
                        else
                            presentEventOrder = presentEventOrder+1;
                            continue
                        end
                end
                
                pw = [ip.dataDir filesep ip.allExperiments{ee} filesep tsDirs{ts} filesep 'eventROIs' filesep cellDirsField{iCellDirs} filesep cellEventsField{iCellEvents}];
                cd (char(pw));
                
                clear roiMetadata
                load('roiMetadata.mat');
                nROIs = nROIs + 1;
                
%               Add metadata to panelMetadata  and labels    
                panelMetadata(nROIs,1) = {char(ip.allExperiments(ee))};                   % EXP000
                panelMetadata(nROIs,2) = {char(tsDirs(ts))};                              % EXP000_TS00
                panelMetadata(nROIs,3) = {char(cellDirs(iCellDirs))};                     % cell_0
                panelMetadata(nROIs,4) = {char(cellEvents(iCellEvents))};                 % cell_0_event_0
                panelMetadata(nROIs,5) = {roiMetadata.xEnd - roiMetadata.xLim + 1};       % xDIM
                panelMetadata(nROIs,6) = {roiMetadata.yEnd - roiMetadata.yLim + 1};       % yDim
                panelMetadata(nROIs,7) = {roiMetadata.nPlanes};                           % nPlanes
                panelMetadata(nROIs,8) = {(roiMetadata.tEnd) - (roiMetadata.tStart) + 1}; % nTimes
                panelMetadata(nROIs,9) = {roiMetadata.tCC};                               % tCC
                panelMetadata(nROIs,12) = {tsDirsField{ts}};                              % tsDirsField
                panelMetadata(nROIs,13) = {cellDirsField{iCellDirs}};                     % cellDirsField
                panelMetadata(nROIs,14) = {cellEventsField{iCellEvents}};                 % cellEventsField
                panelMetadata(nROIs,15) = {roiMetadata.tStart};                           % tStart
                panelMetadata(nROIs,16) = {roiMetadata.tEnd};                             % tEnd
                panelMetadata(nROIs,17) = {roiMetadata.allEventFramesIndex};              % allEventFramesIndex
                panelMetadata(nROIs,18) = {roiMetadata.prevExcludedEvents};               % previously ExcludedEvents in cell
                
                panelMetadataLabels = { ...
                'EXP000'         ; ...
                'EXP000_TS00'    ; ...
                'cell_0'         ; ...
                'cell_0_event_0' ; ...
                'xDim'           ; ...
                'yDim'           ; ...
                'nPlanes'        ; ...
                'nTimes'         ; ...
                'tCC'            ; ...
                'xLim'           ; ...
                'yLim'           };
                panelMetadataLabels = panelMetadataLabels'; %#ok<NASGU>
                
                if roiMetadata.nPlanes > dMaxPlanes
                   dMaxPlanes = roiMetadata.nPlanes;
                end
                
                if (roiMetadata.tEnd - roiMetadata.tStart + 1) > dMaxTimes
                   dMaxTimes = (roiMetadata.tEnd - roiMetadata.tStart + 1);
                end
                
                if roiMetadata.tCC > dMax_tCC
                   dMax_tCC = roiMetadata.tCC;
                end
                cd ../
            end
            clear cellEvents cellEventsField
            cd ../
        end
        clear cellDirs cellDirsField
        cd ../..
    end
    clear tsDirs tsDirsField
    cd ../
end

if isempty(panelMetadata) == 1
    fprintf('\nAll events have been validated.\n')
    return
end

% Start creating panel
panelLayout = 0;
startROI = 1;
panelrow = 1;
while panelfilled == 0
    
    % Determine which events to put in current panel row    
    override = 0;
    iROI = startROI;
    cumSizeX = ip.spacersx;
        
    while override == 0
          if iROI > nROIs
             override = 1; %#ok<NASGU>
             break
          end

          for nROIs = startROI:size(panelMetadata,1)

              xDim = cell2mat(panelMetadata(nROIs,5));
              cumSizeX = cumSizeX + xDim + ip.spacersx;

              if cumSizeX > ip.frameSizeX
                 override = 2;
                 break
              else
                 lastPanelLayout = panelLayout; %#ok<NASGU>
                 panelLayout(panelrow, iROI-startROI+1) = iROI;
                 iROI = iROI + 1;
                 fitSize = cumSizeX;
              end
          end
    end
         
    % Rewind parameters one step
    lastROI = iROI-1;

    blanks(panelrow) = ip.frameSizeX - fitSize + (lastROI-startROI+2)*ip.spacersx;  %#ok<*AGROW>
    spacerxmod(panelrow) = floor(blanks(panelrow)/(lastROI-startROI+2));     
    spacerymod(panelrow) = cumSizeX;                                          %#ok<NASGU>

    % Update lower bound of row
    yDim_row = cell2mat(panelMetadata(startROI:lastROI,6));
    dMax_yDim_row(panelrow) = max(yDim_row);
    cumSizeX = cumSizeX + dMax_yDim_row + ip.spacersy; %#ok<NASGU>

    % Start next panel row of crops remain
    panelrow = panelrow+1;
    if lastROI == size(panelMetadata,1)
       panelfilled = 1; 
    else
       startROI = lastROI+1;
    end
end

% Calculate bounds of all ROIs
[A, B] = size(panelLayout);
yLimCum = 0;
for iA = 1:A
    yLimCum = yLimCum + ip.spacersy;
    xLimCum = 0;
    for iB = 1:B
        idROI = panelLayout(iA, iB);
        if idROI ~= 0
            xLimCum = xLimCum + spacerxmod(iA);
            xDim = panelMetadata(idROI, 5);
            panelMetadata(idROI,10) = {xLimCum};       % xLim of roi
            panelMetadata(idROI,11) = {yLimCum};       % yLim of roi
            xLimCum = xLimCum + cell2mat(xDim);
        end
    end
    yLimCum = yLimCum + dMax_yDim_row(iA);
end

% Add ROIs to current panel row
switch ip.mode
    case 'fullUntracked'
        panel = zeros(ip.frameSizeY, ip.frameSizeX, dMaxPlanes, dMaxTimes, nChannels);
    case 'croppedTracked'
        panel = zeros(ip.frameSizeY, ip.frameSizeX, ip.extPlanes*2+1, dMaxTimes, nChannels+1);
    case 'mipUntracked'
        panel = zeros(ip.frameSizeY, ip.frameSizeX, 1, dMaxTimes, nChannels);
    case 'mipTracked'
        panel = zeros(ip.frameSizeY, ip.frameSizeX, 1, dMaxTimes, nChannels+1);
end

for iROI = 1:size(panelMetadata,1)
    fprintf('.');
    timesAdded = 0;
    xLim = cell2mat(panelMetadata(iROI,10));
    yLim = cell2mat(panelMetadata(iROI,11));
    xDim = cell2mat(panelMetadata(iROI,5));
    yDim = cell2mat(panelMetadata(iROI,6));
    nPlanes = cell2mat(panelMetadata(iROI,7));
    nTimes = cell2mat(panelMetadata(iROI,8));
    exp = cell2mat(panelMetadata(iROI,1));
    ts = cell2mat(panelMetadata(iROI,2));
    CELL = cell2mat(panelMetadata(iROI,3)); %#ok<NASGU>
    event = cell2mat(panelMetadata(iROI,4)); %#ok<NASGU>
    tsField = cell2mat(panelMetadata(iROI,12)); %#ok<NASGU>
    cellField = cell2mat(panelMetadata(iROI,13));
    eventField = cell2mat(panelMetadata(iROI,14));
    firstFrame = cell2mat(panelMetadata(iROI,15));

    pw = [ip.dataDir filesep exp filesep ts filesep 'eventROIs' filesep cellField filesep eventField];
    cd (char(pw));
    load('eventROI.mat')
    load('roiMetadata.mat')
    
    allChannels = ip.allChannels.(char(panelMetadata(iROI,1)));
    nChannels = length(allChannels);
    cd (char(ip.fdp))
    
    switch ip.mode
        case {'croppedTracked', 'mipTracked'}
            pw = [ip.dataDir filesep exp filesep ts filesep 'cargoTrackingdata' filesep cellField filesep eventField];
            try
                cd (char(pw));
            catch me
                continue
            end
                
            load('tracks.mat')
            load('trackingdata.mat')
            tCC = trackingdata.t0;
            
        case {'fullUntracked', 'mipUntracked'}
            tCC = cell2mat(panelMetadata(iROI,9));
    end
    
    switch ip.mode
        case 'fullUntracked'
            for iTime = 1:dMaxTimes
                
                if dMax_tCC-iTime < cell2mat(panelMetadata(iROI,9)) && ... %tCC
                        timesAdded < size(eventROI,4) %#ok<*NODEF>
                    
                    timesAdded = timesAdded + 1;
                    
                    panel(yLim:yLim+yDim-1,xLim:xLim+xDim-1,1:nPlanes,iTime,:) = ...
                        eventROI(:,:,:,timesAdded,:);
                else
                    panel(yLim:yLim+yDim-1,xLim:xLim+xDim-1,1:nPlanes,iTime,:) = 2^16;
                    
                end
                
                % Add text
                for iPlanes = 1:nPlanes
                    normTime = timesAdded - tCC;
                    textbox = zeros(75,170);
                    hh = size(textbox,1);
                    positions = [0 hh-50; 0 hh-35; 0 hh-20; 0 hh-5];
                    textTime = sprintf('Time:  %d', normTime);
                    if timesAdded == 0 || iTime > dMax_tCC + nTimes - tCC
                        normTime = '*';
                        textTime = sprintf('Time:  %s', normTime);
                    end
                    textPlane = sprintf('Plane: %d', iPlanes);
                    text = {char(cell2mat(panelMetadata(iROI,2))), char(cell2mat(panelMetadata(iROI,4))), char(textTime), char(textPlane)};
                    textbox = insertText(textbox,positions,text, 'AnchorPoint', 'LeftBottom', 'BoxColor', 'black', 'TextColor', 'white','BoxOpacity', 0);
                    textbox = im2bw(textbox);
                    textbox = textbox*2^16;
                    textbox(:,:,1,1,2)=textbox;
                    panel(yLim-size(textbox,1):yLim-1, xLim:xLim+size(textbox,2)-1, iPlanes, iTime, :) = textbox;
                end
            end
            
        case 'croppedTracked'
            for iTime = 1:dMaxTimes
                if dMax_tCC-iTime < trackingdata.t0 && ... %tCC
                        timesAdded < size(eventROI,4)
                    
                    timesAdded = timesAdded + 1;
                    firstTrackedTime = trackingdata.coordinates.data(1,4);
                    if timesAdded < firstTrackedTime
                        centerPlane = trackingdata.coordinates.data(1,3);
                    elseif timesAdded >= firstTrackedTime && ...
                            timesAdded-firstTrackedTime+1 <= size(trackingdata.coordinates.data, 1)
                        centerPlane = trackingdata.coordinates.data(timesAdded-firstTrackedTime+1,3);
                    elseif timesAdded-firstTrackedTime+1 > size(trackingdata.coordinates.data, 1)
                        centerPlane = trackingdata.coordinates.data(size(trackingdata.coordinates.data, 1),3);
                    end
                    
                    panelPlane = 0;
                    for addPlanes = centerPlane-ip.extPlanes:centerPlane+ip.extPlanes
                        panelPlane = panelPlane + 1;
                        if addPlanes < 1 || addPlanes > nPlanes
                            panel(yLim:yLim+yDim-1,xLim:xLim+xDim-1,panelPlane,iTime,1:nChannels) = 0;
                            
                        else
                            panel(yLim:yLim+yDim-1,xLim:xLim+xDim-1,panelPlane,iTime,1:nChannels) = ...
                                eventROI(:,:,addPlanes,timesAdded,:);
                            
                            panel(yLim:yLim+yDim-1,xLim:xLim+xDim-1,ip.extPlanes+1,iTime,nChannels+1) = ...
                                tracks(:,:,centerPlane, timesAdded);
                        end
                    end
                    
                else
                    panel(yLim:yLim+yDim-1,xLim:xLim+xDim-1,1:ip.extPlanes*2+1,iTime,1:nChannels) = 2^16;
                    centerPlane = 0;
                end
                
                % Add text
                for iPlanes = 1:ip.extPlanes*2+1
                    planeNumbers = num2cell(centerPlane-ip.extPlanes:centerPlane+ip.extPlanes);
                    
                    normTime = timesAdded - tCC;
                    textbox = zeros(75,170);
                    hh = size(textbox,1);
                    positions = [0 hh-50; 0 hh-35; 0 hh-20; 0 hh-5];
                    textTime = sprintf('Time:  %d', normTime);
                    if timesAdded == 0 || iTime > dMax_tCC + nTimes - tCC
                        normTime = '*';
                        textTime = sprintf('Time:  %s', normTime);
                    end
                    if planeNumbers{iPlanes} < 1
                        textPlane = 'Plane: MIN';
                    elseif planeNumbers{iPlanes} > nPlanes
                        textPlane = 'Plane: MAX';
                    else
                        textPlane = sprintf('Plane: %d', cell2mat(planeNumbers(iPlanes)));
                    end
                    
                    text = {char(cell2mat(panelMetadata(iROI,2))), char(cell2mat(panelMetadata(iROI,4))), char(textTime), char(textPlane)};
                    textbox = insertText(textbox,positions,text, 'AnchorPoint', 'LeftBottom', 'BoxColor', 'black', 'TextColor', 'white','BoxOpacity', 0);
                    textbox = im2bw(textbox); %#ok<*IM2BW>
                    textbox = textbox*2^16;
                    textbox(:,:,1,1,nChannels)=textbox;
                    panel(yLim-size(textbox,1):yLim-1, xLim:xLim+size(textbox,2)-1, iPlanes, iTime,1:nChannels) = textbox;
                end
            end
            
        case 'mipTracked'
            for iTime = 1:dMaxTimes
                if dMax_tCC-iTime < trackingdata.t0 && ... %tCC
                        timesAdded < size(eventROI,4)
                    
                    timesAdded = timesAdded + 1;
                    firstTrackedTime = trackingdata.coordinates.data(1,4);
                    if timesAdded < firstTrackedTime
                        centerPlane = trackingdata.coordinates.data(1,3);
                    elseif timesAdded >= firstTrackedTime && ...
                            timesAdded-firstTrackedTime+1 <= size(trackingdata.coordinates.data, 1)
                        centerPlane = trackingdata.coordinates.data(timesAdded-firstTrackedTime+1,3);
                    elseif timesAdded-firstTrackedTime+1 > size(trackingdata.coordinates.data, 1)
                        centerPlane = trackingdata.coordinates.data(size(trackingdata.coordinates.data, 1),3);
                    end
                    
                    panelPlane = 0;
                    minPlane = centerPlane-ip.extPlanes;
                    maxPlane = centerPlane+ip.extPlanes;
                    if minPlane < 1
                        subROI = zeros(yDim,xDim,(ip.extPlanes*2 + minPlane), 1, nChannels);
                    elseif maxPlane > nPlanes
                        subROI = zeros(yDim, xDim,((ip.extPlanes*2) - (maxPlane-nPlanes)), 1, nChannels);
                    else
                        subROI = zeros(yDim, xDim,(ip.extPlanes*2+1), 1, nChannels);
                    end
                    
                    iSubROI = 1;
                    for addPlanes = centerPlane-ip.extPlanes:centerPlane+ip.extPlanes
                        panelPlane = panelPlane + 1;
                        if addPlanes < 1 || addPlanes > nPlanes
                            continue    
                        else
                            subROI(:,:,iSubROI,1,:) = eventROI(:,:,addPlanes,timesAdded,:);
                            iSubROI = iSubROI+1;
                        end
                    end     
                    
                    MIPsubROI = max(subROI,[],3);
                    
                    panel(yLim:yLim+yDim-1,xLim:xLim+xDim-1, 1,iTime,1:nChannels) = ...
                        MIPsubROI;
                    
                    panel(yLim:yLim+yDim-1,xLim:xLim+xDim-1,1,iTime,nChannels+1) = ...
                        tracks(:,:,centerPlane, timesAdded)*10^3;
                        
                else
                    panel(yLim:yLim+yDim-1,xLim:xLim+xDim-1,1,iTime,1:nChannels) = 2^16;
                    centerPlane = 0;
                end
                
                % Add text    
                normTime = timesAdded - tCC;
                textbox = zeros(75,170);
                hh = size(textbox,1);
                positions = [0 hh-50; 0 hh-35; 0 hh-20; 0 hh-5];
                textTime = sprintf('Time:  %d', normTime);
                if timesAdded == 0 || iTime > dMax_tCC + nTimes - tCC
                    normTime = '*';
                    textTime = sprintf('Time:  %s', normTime);
                end

                textPlane = 'MIP';

                text = {char(cell2mat(panelMetadata(iROI,2))), char(cell2mat(panelMetadata(iROI,4))), char(textTime), char(textPlane)};
                textbox = insertText(textbox,positions,text, 'AnchorPoint', 'LeftBottom', 'BoxColor', 'black', 'TextColor', 'white','BoxOpacity', 0);
                textbox = im2bw(textbox);
                textbox = textbox*2^16;

                if nChannels == 2
                   textbox(:,:,1,1,2)=textbox;
                end
                if nChannels == 3
                    expTextbox = zeros(75,170,1,1,3);
                    for ii = 1:3
                        expTextbox(:,:,:,:,ii) = textbox;
                    end
                   textbox = expTextbox;
                end                
                
                panel(yLim-size(textbox,1):yLim-1, xLim:xLim+size(textbox,2)-1, 1, iTime,1:nChannels) = textbox;
            end
            
            
        case 'mipUntracked'
            firstFlag = 0;
            for iTime = 1:dMaxTimes
                
                if dMax_tCC-iTime < cell2mat(panelMetadata(iROI,9)) && ... %tCC
                        timesAdded < size(eventROI,4)
                    timesAdded = timesAdded + 1;
                    
                    if firstFlag == 0 || timesAdded == size(eventROI,4)
                        
                        panel(yLim:yLim+yDim-1,xLim:xLim+xDim-1, 1,iTime,:) = ...
                        eventROI(:,:,1,timesAdded,:);
                        firstFlag = 1;
                    else
                    

                        panel(yLim:yLim+yDim-1,xLim:xLim+xDim-1, 1,iTime,:) = ...
                            mean(eventROI(:,:,1,timesAdded-1:timesAdded+1,:),4);
                                        
                    end
                else
                    panel(yLim:yLim+yDim-1,xLim:xLim+xDim-1, 1,iTime,:) = 2^16;
                    
                end
                
                % Add text
                normTime = timesAdded - tCC;
                frameNumber = firstFrame + timesAdded;
                textbox = zeros(105,170);
                hh = size(textbox,1);
                positions = [0 hh-80; 0 hh-65; 0 hh-50; 0 hh-35; 0 hh-20; 0 hh-5];
                textTime = sprintf('Time:  %d', normTime);
                textFrame = sprintf('Frame:  %d', frameNumber);
                if timesAdded == 0 || iTime > dMax_tCC + nTimes - tCC
                    normTime = '*';
                    textTime = sprintf('Time:  %s', normTime);
                    frameNumber = '*';
                    textFrame = sprintf('Frame:  %s', frameNumber);
                end
                textPlane = 'MIP';
                text = {char(cell2mat(panelMetadata(iROI,2))), char(cell2mat(panelMetadata(iROI,4))), char(textTime), char(textFrame),...
                    char(sprintf('allEventFrames# %d', panelMetadata{iROI,17})), char(sprintf('cellEvent# %d', panelMetadata{iROI,18}+1))};
                textbox = insertText(textbox,positions,text, 'AnchorPoint', 'LeftBottom', 'BoxColor', 'black', 'TextColor', 'white','BoxOpacity', 0);
                textbox = im2bw(textbox);
                textbox = textbox*2^16;
                if nChannels == 2
                   textbox(:,:,1,1,2)=textbox;
                end
                if nChannels == 3
                    expTextbox = zeros(105,170,1,1,3);
                    for ii = 1:3
                        expTextbox(:,:,:,:,ii) = textbox;
                    end
                   textbox = expTextbox;
                end
                panel(yLim-size(textbox,1):yLim-1, xLim:xLim+size(textbox,2)-1, 1, iTime, :) = textbox;
            end
            
   end
end

% Adjust ip.frameSizeY to the size of ROIs (last row)
if yLimCum+ip.spacersy > ip.frameSizeY
   extend = (yLimCum+ip.spacersy)-(size(panel,1)-1);
   panel(size(panel,1):size(panel,1)+extend,:,:,:,:) = 0;
elseif (yLimCum+ip.spacersy)-ip.frameSizeY < 0
       panel((yLimCum+ip.spacersy):size(panel,1),:,:,:,:) = [];
end

cd (ip.fdp)
exportEventPanel(ip,panel,panelMetadata);

fprintf('\nEvent panel has been exported.\n');

end

