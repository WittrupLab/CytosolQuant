function [] = secondaryProcessing(ip)
%% Secondary processing operations
fprintf('Secondary processing...\n');
allExperiments = ip.allExperiments;
nExperiments = length(allExperiments);

% Create allROI metadata struct
cd (char(ip.dataDir))

for ee = 1:nExperiments
    pw = [ip.dataDir filesep allExperiments{ee}];
    cd (pw);
    
    dirData = dir; % Get the data for the current directory
    dirData(strncmp({dirData.name}, '.', 1)) = []; % Remove false files
    dirData(strncmp({dirData.name}, '~', 1)) = [];
    dirIndex = [dirData.isdir];  % Find the index for directories
    tsDirs = {dirData(dirIndex).name}';
    tsDirsField = tsDirs;
    pos1Dir = tsDirs{1};
    
    cd (char(pos1Dir));
    cdir = pwd;
    cd (ip.fdp);
    
    % use function for calculation
    illumReference = calcIllumReference(ip,cdir);
    path = [ip.dataDir filesep (allExperiments{ee})];
    cd (path)
    
    clear nonEventCells controlCells
    nonEventCells = struct;
    controlCells = struct;
    
    index_1 = 1;
    index_2 = 1;

    for TS = 1:length(tsDirs)
        cd (tsDirs{TS});

        % Load cellDataStructure file
        cd 'matlabOutput'
        clear dir
        dirData = dir; % Get the data for the current directory
        dirData(strncmp({dirData.name}, '.', 1)) = []; % Remove false files
        dirData(strncmp({dirData.name}, '~', 1)) = [];
        allFiles = {dirData.name}';
        
        for uu = 1:length(allFiles)
            if contains (allFiles(uu),'cellDataStructure') == 1
                clear cellDataStructure
                load(char(allFiles(uu)));
            end
        end
        
        experiment = allExperiments{ee};
        position = tsDirsField{TS};
        
        allCells = fieldnames(cellDataStructure.(allExperiments{ee}).(tsDirsField{TS}));
        for cc = 1:length(allCells)
            cell = allCells{cc};
            
            try dd = cellDataStructure.(experiment).(position).(cell).disqualifyCell;
            catch me
                cellDataStructure.(experiment).(position).(cell).disqualifyCell = 'no';
            end
            
            switch cellDataStructure.(experiment).(position).(cell).disqualifyCell
                case 'yes'
                    continue
            end
            
            cellEvents = fieldnames(cellDataStructure.(experiment).(position).(cell).allEvents);
            % Default nonEventCell = yes, deny if event present
            nonEventCell = 'yes';
            
            switch cellDataStructure.(experiment).(position).(cell).eGFPExpressionLevel
                case 'subthreshold'
                    continue
            end
            
            if isnan(cellDataStructure.(experiment).(position).(cell).eGFPExpressionLevel)
                continue
            end
            
            if cellDataStructure.(experiment).(position).(cell).cellTrackingStarted > 3
                continue
            end
            
            if cellDataStructure.(experiment).(position).(cell).trackDuration < 50
                continue
            end
            
            switch cellDataStructure.(experiment).(position).(cell).controlPosition
                case 'yes'
                    maskTraceAt = Inf;
                    cellDataStructure.(experiment).(position).(cell).controlCell = 'yes';
                    
                    traces = {'rawTraces', 'mfTraces', 'bcTraces', 'mfBcTraces'};
                    channels = {'eGFP', 'siRNA', 'siRNAF'};
                    
                    apoptosisFrame = cellDataStructure.(experiment).(position).(cell).apoptosis.frame;
                    if isnan(apoptosisFrame) ==  0 %ei detected
                        % dont collect if cell dies
                        % within 30 frames from
                        % event. If apoptosis is
                        % present mask 20 frames
                        % prior
                        if apoptosisFrame < 30
                            break
                        else
                            maskTraceAt = apoptosisFrame - 20;
                        end
                    end
                    
                    for iTraces = 1:numel(traces)
                        for iChannels = 1:numel(channels)
                            
                            controlCells.data.(traces{iTraces}).(channels{iChannels})(1:1000,index_1) = NaN;
                            data = cellDataStructure.(experiment).(position).(cell).(traces{iTraces}).(channels{iChannels});
                            controlCells.data.(traces{iTraces}).(channels{iChannels})(1:length(data),index_1) = data;
                            % MASK DYING CELLS
                            maskedData = data;
                            if maskTraceAt < Inf
                                maskedData(maskTraceAt:end) = NaN;
                            end
                            nonEventCells.data.(traces{iTraces}).(channels{iChannels})(1:length(data),index_2) = maskedData;
                            
                        end
                    end
                    controlCells.ID(index_1,1) = {sprintf('Column_%d',index_1)};
                    controlCells.ID(index_1,2) = {experiment};
                    controlCells.ID(index_1,3) = {position};
                    controlCells.ID(index_1,4) = {cell};
                    index_1 = index_1 + 1;
                    continue
            end
            
            for iCellEvents = 1:length(cellEvents)
                validationStatus = cellDataStructure.(experiment).(position).(cell).allEvents.(cellEvents{iCellEvents}).validationStatus;
                switch validationStatus
                    case 'not validated'
                        break
                        
                    case 'false event'
                        continue
                        
                    case 'valid event'
                        nonEventCell = 'no';
                        eventFrame = cellDataStructure.(experiment).(position).(cell).allEvents.(cellEvents{iCellEvents}).frame;
                        xCoord = floor(cellDataStructure.(experiment).(position).(cell).trackingCoordinates.yCoordinate(eventFrame,1));
                        yCoord = floor(cellDataStructure.(experiment).(position).(cell).trackingCoordinates.xCoordinate(eventFrame,1));
                        illumReferenceVar = illumReference(xCoord,yCoord);
                        cellDataStructure.(experiment).(position).(cell).allEvents.(cellEvents{iCellEvents}).illumReferenceVar = illumReferenceVar;
                end
            end % iCellEvents
            
            cellDataStructure.(experiment).(position).(cell).nonEventCell = nonEventCell;
            switch nonEventCell
                case 'yes'
                    maskTraceAt = Inf;
                    traces = {'rawTraces', 'mfTraces', 'bcTraces', 'mfBcTraces'};
                    channels = {'eGFP', 'siRNA', 'siRNAF'};
                    
                    apoptosisFrame = cellDataStructure.(experiment).(position).(cell).apoptosis.frame;
                    if isnan(apoptosisFrame) ==  0 %ei detected
                        % dont collect if cell dies
                        % within 30 frames from
                        % event. If apoptosis is
                        % present mask 20 frames
                        % prior
                        if apoptosisFrame < 30
                            break
                        else
                            maskTraceAt = apoptosisFrame - 20;
                        end
                    end
                    
                    for iTraces = 1:numel(traces)
                        for iChannels = 1:numel(channels)

                            nonEventCells.data.(traces{iTraces}).(channels{iChannels})(1:1000,index_2) = NaN;
                            data = cellDataStructure.(experiment).(position).(cell).(traces{iTraces}).(channels{iChannels});
                            % MASK DYING CELLS
                            maskedData = data;
                            if maskTraceAt < Inf
                                maskedData(maskTraceAt:end) = NaN;
                            end
                            nonEventCells.data.(traces{iTraces}).(channels{iChannels})(1:length(data),index_2) = maskedData;
                            
                        end
                    end
                    nonEventCells.ID(index_2,1) = {sprintf('Column_%d',index_2)};
                    nonEventCells.ID(index_2,2) = {experiment};
                    nonEventCells.ID(index_2,3) = {position};
                    nonEventCells.ID(index_2,4) = {cell};
                    index_2 = index_2 + 1;   
            end
        end
        filename = char(sprintf('%s%s', position,'_cellDataStructure'));
        save(char(filename),'cellDataStructure')
        cd ../.. 
        
    end %tsDirs
    
    if numel(fieldnames(controlCells)) == 0
        controlCells = 'No existing control positions'; %#ok<*NASGU>
    end
    
    for TS = 1:length(tsDirs)
        cd (tsDirs{TS});
        position = tsDirsField{TS};
        
        %  Save nonEventCells and controlCells data files
        cd 'matlabOutput'

        filename = char(sprintf('%s%s', position,'_controlCells'));
        save(char(filename),'controlCells')
        filename = char(sprintf('%s%s', position,'_nonEventCells'));
        save(char(filename),'nonEventCells')
        cd ../..
    end
end

fprintf('\nSecondary processing completed.\n');
end

