function [eventCells_primary] = collectData(ip)
%% Compose time series panel with all events and tracks of the capture
fprintf('\nCollecting data...\n')
allExperiments = ip.allExperiments;
nExperiments = length(allExperiments);
traces = {'rawTraces', 'mfTraces', 'mfBcTraces' 'bcTraces'};
channels = {'eGFP', 'siRNA', 'siRNAF'};

%nNonEventCells
nControlCells = struct;
EventCells = struct;
eventCells_primary = struct;

indexMitosis_primary = 1;
indexMitosis_controlCells = 1;
indexMitosis_nonEventCells = 1;
index_1 = 1;
index_2 = 1;
index_3 = 1;
index_4 = 1;
index_5 = 1;
index_6 = 1;
index_7 = 1;
index_8 = 1;
index_9 = 1;
index_10 = 1;
index_11 = 1;
index_12 = 1;

% Create allROI metadata struct
cd (char(ip.dataDir))

for ee = 1:nExperiments
    cd (char(allExperiments(ee)))
    if ip.processAll_TS == 1
        dirData = dir; % Get the data for the current directory
        dirData(strncmp({dirData.name}, '.', 1)) = []; % Remove false files
        dirData(strncmp({dirData.name}, '~', 1)) = [];
        dirIndex = [dirData.isdir];  % Find the index for directories
        tsDirs = {dirData(dirIndex).name}';
        tsDirsField = tsDirs;
    elseif ip.processAll_TS == 0
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
    
    for TS = 1:length(tsDirs)
        cd (char(tsDirs(TS)));
        
        % load detectedReleaseEvents file
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
            elseif contains (allFiles(uu),'controlCells') == 1
                clear controlCells
                load(char(allFiles(uu)));
            elseif contains (allFiles(uu),'nonEventCells') == 1
                clear nonEventCells
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
                    
                    if ip.collect.controlCells == 1
                        maskTraceAt = Inf;
                        cellTrackingStarted = cellDataStructure.(experiment).(position).(cell).cellTrackingStarted;
                        cellTrackingEnded = cellDataStructure.(experiment).(position).(cell).cellTrackingEnded;
                        
                        apoptosisFrame = cellDataStructure.(experiment).(position).(cell).apoptosis.frame;
                        if isnan(apoptosisFrame) ==  0 %ie detected
                            % dont collect if cell dies
                            % within 30 frames from
                            % event. If apoptosis is
                            % present mask 20 frames
                            % prior
                            if apoptosisFrame - cellTrackingStarted < 30
                                continue
                            else
                                maskTraceAt = apoptosisFrame - 20;
                            end
                        end
                        
                        for iTraces = 1:numel(traces)
                            for iChannels = 1:numel(channels)
                                nControlCells.data.(traces{iTraces}).(channels{iChannels})(1:1000,index_1) = NaN;
                                nControlCells.maskedData.(traces{iTraces}).(channels{iChannels})(1:1000,index_1) = NaN;
                                data = cellDataStructure.(experiment).(position).(cell).(traces{iTraces}).(channels{iChannels});
                                
                                if iChannels == 1 %GFP
                                    
                                    % norm to first timepoint
                                    data(cellTrackingStarted:cellTrackingEnded) = ...
                                        data(cellTrackingStarted:cellTrackingEnded) / data(cellTrackingStarted);
                                    
                                else % siRNA
                                    % skip siRNA illumination correction, no events to use
                                end
                                
                                nControlCells.data.(traces{iTraces}).(channels{iChannels})(1:length(data),index_1) = data;
                                
                                % MASK DYING CELLS
                                maskedData = data;
                                if maskTraceAt < Inf
                                    maskedData(maskTraceAt:end) = NaN;
                                end
                                nControlCells.maskedData.(traces{iTraces}).(channels{iChannels})(1:length(data),index_1) = maskedData;
                                
                            end
                        end
                        
                        nonTSdata = cellDataStructure.(experiment).(position).(cell).rawTraces.eGFP;
                        for qq = 1:length(nonTSdata)
                            if isnan(nonTSdata(qq)) == 0
                                nonTSdata_first = qq;
                                break
                            end
                        end
                        
                        mitosisFrame = cellDataStructure.(experiment).(position).(cell).mitosis.first.frame;
                        allMitosis_controlCells(indexMitosis_controlCells) = mitosisFrame;
                        indexMitosis_controlCells = indexMitosis_controlCells + 1;
                        
                        nControlCells.ID(index_1,1) = {sprintf('Column_%d',index_1)};
                        nControlCells.ID(index_1,2) = {experiment};
                        nControlCells.ID(index_1,3) = {position};
                        nControlCells.ID(index_1,4) = {cell};
                        index_1 = index_1 + 1;
                    end
                    
                    continue
                    
            end
            
            try dd = cellDataStructure.(experiment).(position).(cell).nonEventCell;
            catch me
                continue
            end
            
            switch cellDataStructure.(experiment).(position).(cell).nonEventCell
                case 'yes'
                    %perform collecting if selected
                    if ip.collect.nonEventCells == 1
                        maskTraceAt = Inf;
                        cellTrackingStarted = cellDataStructure.(experiment).(position).(cell).cellTrackingStarted;
                        cellTrackingEnded = cellDataStructure.(experiment).(position).(cell).cellTrackingEnded;
                        
                        apoptosisFrame = cellDataStructure.(experiment).(position).(cell).apoptosis.frame;
                        if isnan(apoptosisFrame) ==  0 %ei detected
                            % dont collect if cell dies
                            % within 30 frames from
                            % event. If apoptosis is
                            % present mask 20 frames
                            % prior
                            if apoptosisFrame - cellTrackingStarted < 30
                                continue
                            else
                                maskTraceAt = apoptosisFrame - 20;
                            end
                        end
                        
                        for iTraces = 1:numel(traces)
                            for iChannels = 1:numel(channels)
                                nNonEventCells.data.(traces{iTraces}).(channels{iChannels})(1:1000,index_2) = NaN;
                                nNonEventCells.maskedData.(traces{iTraces}).(channels{iChannels})(1:1000,index_2) = NaN;
                                data = cellDataStructure.(experiment).(position).(cell).(traces{iTraces}).(channels{iChannels});
                                
                                if iChannels == 1 %GFP
                                    data(cellTrackingStarted:cellTrackingEnded) = ...
                                        data(cellTrackingStarted:cellTrackingEnded) / data(cellTrackingStarted); %norm to first timepoint
                                else %siRNA
                                    %                                     skip siRNA illumination correction, no events to use
                                end
                                
                                nNonEventCells.data.(traces{iTraces}).(channels{iChannels})(1:length(data),index_2) = data;
                                
                                % MASK DYING CELLS
                                maskedData = data;
                                if maskTraceAt < Inf
                                    maskedData(maskTraceAt:end) = NaN;
                                end
                                nNonEventCells.maskedData.(traces{iTraces}).(channels{iChannels})(1:length(data),index_2) = maskedData;
                                
                            end
                        end
                        
                        mitosisFrame = cellDataStructure.(experiment).(position).(cell).mitosis.first.frame;
                        allMitosis_nonEventCells(indexMitosis_nonEventCells) = mitosisFrame;
                        indexMitosis_nonEventCells = indexMitosis_nonEventCells + 1;
                        
                        nNonEventCells.ID(index_2,1) = {sprintf('Column_%d',index_2)};
                        nNonEventCells.ID(index_2,2) = {experiment};
                        nNonEventCells.ID(index_2,3) = {position};
                        nNonEventCells.ID(index_2,4) = {cell};
                        index_2 = index_2 + 1;
                    end
                    
                    if ip.collect.mitosisCells_nonEvenCells_timeShifted == 1 %mitosis at '0'
                        cellMitosis = fieldnames(cellDataStructure.(experiment).(position).(cell).mitosis);
                        for iMitosis = 1:numel(cellMitosis)
                            mitosisStatus = cellDataStructure.(experiment).(position).(cell).mitosis.(cellMitosis{iMitosis}).status;
                            switch mitosisStatus
                                case 'detected'
                                    for iTraces = 1:numel(traces)
                                        for iChannels = 1:numel(channels)
                                            mitosis_nonEventCells.(cellMitosis{iMitosis}).data.(traces{iTraces}).(channels{iChannels})(1:1000,index_3) = NaN;
                                            data = cellDataStructure.(experiment).(position).(cell).tsTraces.mitosis.(cellMitosis{iMitosis}).(traces{iTraces}).(channels{iChannels});
                                            mitosis_nonEventCells.(cellMitosis{iMitosis}).data.(traces{iTraces}).(channels{iChannels})(1:length(data),index_3) = data;
                                        end
                                    end
                                    mitosis_nonEventCells.(cellMitosis{iMitosis}).ID(index_3,1) ={sprintf('Column_%d',index_3)};
                                    mitosis_nonEventCells.(cellMitosis{iMitosis}).ID(index_3,2) = {experiment};
                                    mitosis_nonEventCells.(cellMitosis{iMitosis}).ID(index_3,3) = {position};
                                    mitosis_nonEventCells.(cellMitosis{iMitosis}).ID(index_3,4) = {cell};
                                    index_3 = index_3 + 1;
                            end
                            
                        end
                    end
                    
                    if ip.collect.apoptoticNonEventCells == 1 % apoptosis at '0'
                        apoptosisStatus = fieldnames(cellDataStructure.(experiment).(position).(cell).apoptosis.status);
                        switch apoptosisStatus
                            case 'detected'
                                for iTraces = 1:numel(traces)
                                    for iChannels = 1:numel(channels)
                                        apoptosis_nonEventCells.data.(traces{iTraces}).(channels{iChannels})(1:1000,index_4) = NaN;
                                        data = cellDataStructure.(experiment).(position).(cell).tsTraces.apoptosis.(traces{iTraces}).(channels{iChannels});
                                        apoptosis_nonEventCells.data.(traces{iTraces}).(channels{iChannels})(1:length(data),index_4) = data;
                                    end
                                end
                                apoptosis_nonEventCells.ID(addEvent,1) = sprintf('Column_%d',index_4);
                                apoptosis_nonEventCells.ID(addEvent,2) = {experiment};
                                apoptosis_nonEventCells.ID(addEvent,3) = {position};
                                apoptosis_nonEventCells.ID(addEvent,4) = {cell};
                                index_4 = index_4 + 1;
                        end
                    end
                    
                case 'no'% eventCells
                    
                    if ip.collect.eventCells_timeShifted_primary == 1
                        presentEventOrder = 1;
                        cellEvents = fieldnames(cellDataStructure.(experiment).(position).(cell).allEvents);
                        for iCellEvents = 1:length(cellEvents)
                            
                            validationStatus = cellDataStructure.(experiment).(position).(cell).allEvents.(cellEvents{iCellEvents}).validationStatus;
                            
                            switch validationStatus
                                case 'not validated'
                                    fprintf('\nCannot collect non-validated event, skipping cell...\n')
                                    break
                                    
                                case 'false event'
                                    continue
                                    
                                case 'valid event'
                                    % check temporal proximity to apoptosis
                                    maskTraceAt = Inf;
                                    if presentEventOrder == 1
                                        eventFrame = cellDataStructure.(experiment).(position).(cell).allEvents.(cellEvents{iCellEvents}).frame;
                                        apoptosisFrame = cellDataStructure.(experiment).(position).(cell).apoptosis.frame;
                                        if isnan(apoptosisFrame) ==  0 %ei detected
                                            % dont collect if cell dies
                                            % within 30 frames from
                                            % event. If apoptosis is
                                            % present mask 20 frames
                                            % prior
                                            if apoptosisFrame - eventFrame < 30
                                                break
                                            else
                                                maskTraceAt = apoptosisFrame - 20;
                                            end
                                        end
                                        
                                        %look for next event, if closer
                                        %than apoptosis, change masking
                                        %timepoint
                                        nextValidEventIndex = 1;
                                        for iCellEvents2 = 1:numel(cellEvents)
                                            switch cellDataStructure.(experiment).(position).(cell).allEvents.(cellEvents{iCellEvents2}).validationStatus
                                                case 'valid event'
                                                    if nextValidEventIndex == 2
                                                        nextEventFrame = cellDataStructure.(experiment).(position).(cell).allEvents.(cellEvents{iCellEvents2}).frame;
                                                        if nextEventFrame < maskTraceAt
                                                            maskTraceAt = nextEventFrame;
                                                        end
                                                    else
                                                        nextValidEventIndex = nextValidEventIndex +1;
                                                        continue
                                                    end
                                            end
                                        end
                                        cellTrackingStarted = cellDataStructure.(experiment).(position).(cell).cellTrackingStarted;
                                        cellTrackingEnded = cellDataStructure.(experiment).(position).(cell).cellTrackingEnded;
                                        for iTraces = 1:numel(traces)
                                            for iChannels = 1:numel(channels)
                                                eventCells_primary.noShift.data.(traces{iTraces}).(channels{iChannels})(1:1000,index_5) = NaN;
                                                data = cellDataStructure.(experiment).(position).(cell).(traces{iTraces}).(channels{iChannels});
                                                if iChannels == 1 %GFP
                                                    %bleachCorrFactor = nanmean(controlCells.data.mfBcTraces.eGFP,2);
                                                    bleachCorrFactor = nanmean(nonEventCells.data.mfBcTraces.eGFP,2);
                                                    bleachCorrFactor = bleachCorrFactor ./ nanmean(bleachCorrFactor(1:5));
                                                    data(cellTrackingStarted:cellTrackingEnded) = ...
                                                        data(cellTrackingStarted:cellTrackingEnded) ./ bleachCorrFactor(cellTrackingStarted:cellTrackingEnded);
                                                    
                                                    
                                                    if ip.collect.normGFP == 1
                                                        data(cellTrackingStarted:cellTrackingEnded) = ...
                                                            data(cellTrackingStarted:cellTrackingEnded) ./ data(eventFrame);% normalize to single trace
                                                    end
                                                    
                                                else %siRNA
                                                    illumReferenceVar = cellDataStructure.(experiment).(position).(cell).allEvents.(cellEvents{iCellEvents}).illumReferenceVar;
                                                    data(cellTrackingStarted:cellTrackingEnded) = ...
                                                        data(cellTrackingStarted:cellTrackingEnded) * 2^16 / illumReferenceVar;
                                                    
                                                end
                                                eventCells_primary.noShift.data.(traces{iTraces}).(channels{iChannels})(1:length(data),index_5) = data;
                                                
                                                eventCells_primary.timeShift.data.(traces{iTraces}).(channels{iChannels})(1:1000,index_5) = NaN;
                                                data = cellDataStructure.(experiment).(position).(cell).tsTraces.events.(cellEvents{iCellEvents}).(traces{iTraces}).(channels{iChannels});
                                                if iChannels == 1 %GFP
                                                    %bleachCorrFactor = nanmean(controlCells.data.mfBcTraces.eGFP,2);
                                                    bleachCorrFactor = nanmean(nonEventCells.data.mfBcTraces.eGFP,2);
                                                    bleachCorrFactor = bleachCorrFactor ./ nanmean(bleachCorrFactor(1:5));

                                                    data(500-(eventFrame-cellTrackingStarted):500+(cellTrackingEnded-eventFrame)) = ...
                                                        data(500-(eventFrame-cellTrackingStarted):500+(cellTrackingEnded-eventFrame)) ./ bleachCorrFactor(cellTrackingStarted:cellTrackingEnded);
                                                    
                                                    if ip.collect.normGFP == 1
                                                        data(500-(eventFrame-cellTrackingStarted):500+(cellTrackingEnded-eventFrame)) = ...
                                                            data(500-(eventFrame-cellTrackingStarted):500+(cellTrackingEnded-eventFrame)) ./ data(500);
                                                    end
                                                    
                                                else %siRNA
                                                    illumReferenceVar = cellDataStructure.(experiment).(position).(cell).allEvents.(cellEvents{iCellEvents}).illumReferenceVar;
                                                    data(500-(eventFrame-cellTrackingStarted):500+(cellTrackingEnded-eventFrame)) = ...
                                                        data(500-(eventFrame-cellTrackingStarted):500+(cellTrackingEnded-eventFrame)) * 2^16 / illumReferenceVar;
                                                end
                                                eventCells_primary.timeShift.data.(traces{iTraces}).(channels{iChannels})(1:length(data),index_5) = data;
                                                
                                                maskedData = data;
                                                if maskTraceAt < Inf
                                                    maskedData(500-(eventFrame-cellTrackingStarted)+maskTraceAt:end) = NaN;
                                                end
                                                eventCells_primary.timeShift.maskedData.(traces{iTraces}).(channels{iChannels})(1:length(data),index_5) = maskedData;
                                            end
                                        end
                                        
                                        nonTSdata = cellDataStructure.(experiment).(position).(cell).rawTraces.eGFP;
                                        for qq = 1:length(nonTSdata)
                                            if isnan(nonTSdata(qq)) == 0
                                                nonTSdata_first = qq;
                                                break
                                            end
                                        end
                                        
                                        TSdata = cellDataStructure.(experiment).(position).(cell).tsTraces.events.first.rawTraces.eGFP;
                                        for qq = 1:length(TSdata)
                                            if isnan(TSdata(qq)) == 0
                                                TSdata_first = qq;
                                                break
                                            end
                                        end
                                        
                                        mitosisFrame = cellDataStructure.(experiment).(position).(cell).mitosis.first.frame;
                                        allMitosis_primary(indexMitosis_primary) = mitosisFrame + (TSdata_first - nonTSdata_first);
                                        indexMitosis_primary = indexMitosis_primary + 1;
                                        
                                        eventCells_primary.timeShift.ID(index_5,1) = {sprintf('Column_%d',index_5)};
                                        eventCells_primary.timeShift.ID(index_5,2) = {experiment};
                                        eventCells_primary.timeShift.ID(index_5,3) = {position};
                                        eventCells_primary.timeShift.ID(index_5,4) = {cell};
                                        
                                        eventCells_primary.noShift.ID(index_5,1) = {sprintf('Column_%d',index_5)};
                                        eventCells_primary.noShift.ID(index_5,2) = {experiment};
                                        eventCells_primary.noShift.ID(index_5,3) = {position};
                                        eventCells_primary.noShift.ID(index_5,4) = {cell};
                                        index_5 = index_5 + 1;
                                        
                                        break% no further events from cell
                                    else
                                        presentEventOrder = presentEventOrder + 1;
                                    end
                            end
                        end % iCellEvents
                    end % ip.collect
                    
                    if ip.collect.eventCells_timeShifted_secondary == 1
                        presentEventOrder = 1;
                        cellEvents = fieldnames(cellDataStructure.(experiment).(position).(cell).allEvents);
                        for iCellEvents = 1:length(cellEvents)
                            
                            validationStatus = cellDataStructure.(experiment).(position).(cell).allEvents.(cellEvents{iCellEvents}).validationStatus;
                            
                            switch validationStatus
                                case 'not validated'
                                    fprintf('\nCannot collect non-validated event, skipping cell...\n')
                                    break
                                    
                                case 'false event'
                                    continue
                                    
                                case 'valid event'
                                    if presentEventOrder == 2
                                        
                                        eventFrame = cellDataStructure.(experiment).(position).(cell).allEvents.(cellEvents{iCellEvents}).frame;
                                        apoptosisFrame = cellDataStructure.(experiment).(position).(cell).apoptosis.frame;
                                        if isnan(apoptosisFrame) ==  0 %ei detected
                                            % dont collect if cell dies
                                            % within 30 frames from
                                            % event. If apoptosis is
                                            % present mask 20 frames
                                            % prior
                                            if apoptosisFrame - eventFrame < 30
                                                break
                                            else
                                                maskTraceAt = apoptosisFrame - 20;
                                            end
                                        end
                                        
                                        %look for next event, if closer
                                        %than apoptosis, change masking
                                        %timepoint
                                        nextValidEventIndex = 1;
                                        for iCellEvents2 = 1:numel(cellEvents)
                                            switch cellDataStructure.(experiment).(position).(cell).allEvents.(cellEvents{iCellEvents2}).validationStatus
                                                case 'valid event'
                                                    if nextValidEventIndex == 3
                                                        nextEventFrame = cellDataStructure.(experiment).(position).(cell).allEvents.(cellEvents{iCellEvents2}).frame;
                                                        if nextEventFrame < maskTraceAt
                                                            maskTraceAt = nextEventFrame;
                                                        end
                                                    else
                                                        nextValidEventIndex = nextValidEventIndex +1;
                                                        continue
                                                    end
                                            end
                                        end
                                        
                                        for iTraces = 1:numel(traces)
                                            for iChannels = 1:numel(channels)
                                                
                                                eventCells_secondary.timeShift.data.(traces{iTraces}).(channels{iChannels})(1:1000,index_6) = NaN;
                                                data = cellDataStructure.(experiment).(position).(cell).tsTraces.events.(cellEvents{iCellEvents}).(traces{iTraces}).(channels{iChannels});
                                                % norm
                                                eventCells_secondary.timeShift.data.(traces{iTraces}).(channels{iChannels})(1:length(data),index_6) = data;
                                                maskedData = data;
                                                maskedData(maskTraceAt:end) = NaN;
                                                eventCells_secondary.timeShift.maskedData.(traces{iTraces}).(channels{iChannels})(1:length(data),index_5) = maskedData;
                                            end
                                        end
                                        eventCells_secondary.timeShift.ID(addEvent,1) = sprintf('Column_%d',index_6);
                                        eventCells_secondary.timeShift.ID(addEvent,2) = {experiment};
                                        eventCells_secondary.timeShift.ID(addEvent,3) = {position};
                                        eventCells_secondary.timeShift.ID(addEvent,4) = {cell};
                                        
                                        index_6 = index_6 + 1;
                                        % look for next valid event or
                                        % apoptosis --> mask remaining trax
                                        break
                                    else
                                        presentEventOrder = presentEventOrder + 1;
                                    end
                            end
                        end % iCellEvents
                    end
                    
                    if ip.collect.eventCells_timeShifted_tertiary == 1
                        presentEventOrder = 1;
                        cellEvents = fieldnames(cellDataStructure.(experiment).(position).(cell).allEvents);
                        for iCellEvents = 1:length(cellEvents)
                            
                            validationStatus = cellDataStructure.(experiment).(position).(cell).allEvents.(cellEvents{iCellEvents}).validationStatus;
                            
                            switch validationStatus
                                case 'not validated'
                                    fprintf('\nCannot collect non-validated event, skipping cell...\n')
                                    break
                                    
                                case 'false event'
                                    continue
                                    
                                case 'valid event'
                                    if presentEventOrder == 3
                                        for iTraces = 1:numel(traces)
                                            for iChannels = 1:numel(channels)
                                                eventCells_tertiary.noShift.data.(traces{iTraces}).(channels{iChannels})(1:1000,index_7) = NaN;
                                                data = cellDataStructure.(experiment).(position).(cell).(traces{iTraces}).(channels{iChannels});
                                                % norm
                                                eventCells_tertiary.noShiftdata.(traces{iTraces}).(channels{iChannels})(1:length(data),index_7) = data;
                                                
                                                eventCells_tertiary.timeShift.data.(traces{iTraces}).(channels{iChannels})(1:1000,index_7) = NaN;
                                                data = cellDataStructure.(experiment).(position).(cell).tsTraces.events.(cellEvents{iCellEvents}).(traces{iTraces}).(channels{iChannels});
                                                % norm
                                                eventCells_tertiary.timeShift.data.(traces{iTraces}).(channels{iChannels})(1:length(data),index_7) = data;
                                            end
                                        end
                                        eventCells_tertiary.timeShift.ID(addEvent,1) = sprintf('Column_%d',index_7);
                                        eventCells_tertiary.timeShift.ID(addEvent,2) = {experiment};
                                        eventCells_tertiary.timeShift.ID(addEvent,3) = {position};
                                        eventCells_tertiary.timeShift.ID(addEvent,4) = {cell};
                                        
                                        eventCells_tertiary.noShift.ID(addEvent,1) = sprintf('Column_%d',index_7);
                                        eventCells_tertiary.noShift.ID(addEvent,2) = {experiment};
                                        eventCells_tertiary.noShift.ID(addEvent,3) = {position};
                                        eventCells_tertiary.noShift.ID(addEvent,4) = {cell};
                                        index_7 = index_7 + 1;
                                        %look for next valid event or
                                        %apoptosis --> mask remaining trax
                                        break
                                    else
                                        presentEventOrder = presentEventOrder + 1;
                                    end
                            end
                        end % iCellEvents
                    end
                    
                    if ip.collect.eventCells_timeShifted_quatary == 1
                        presentEventOrder = 1;
                        cellEvents = fieldnames(cellDataStructure.(experiment).(position).(cell).allEvents);
                        for iCellEvents = 1:length(cellEvents)
                            
                            validationStatus = cellDataStructure.(experiment).(position).(cell).allEvents.(cellEvents{iCellEvents}).validationStatus;
                            
                            switch validationStatus
                                case 'not validated'
                                    fprintf('\nCannot collect non-validated event, skipping cell...\n')
                                    break
                                    
                                case 'false event'
                                    continue
                                    
                                case 'valid event'
                                    if presentEventOrder == 3
                                        for iTraces = 1:numel(traces)
                                            for iChannels = 1:numel(channels)
                                                eventCells_quatary.noShift.data.(traces{iTraces}).(channels{iChannels})(1:1000,index_8) = NaN;
                                                data = cellDataStructure.(experiment).(position).(cell).(traces{iTraces}).(channels{iChannels});
                                                eventCells_quatary.noShiftdata.(traces{iTraces}).(channels{iChannels})(1:length(data),index_8) = data;
                                                
                                                eventCells_quatary.timeShift.data.(traces{iTraces}).(channels{iChannels})(1:1000,index_8) = NaN;
                                                data = cellDataStructure.(experiment).(position).(cell).tsTraces.events.(cellEvents{iCellEvents}).(traces{iTraces}).(channels{iChannels});
                                                eventCells_quatary.timeShift.data.(traces{iTraces}).(channels{iChannels})(1:length(data),index_8) = data;
                                            end
                                        end
                                        eventCells_quatary.timeShift.ID(addEvent,1) = sprintf('Column_%d',index_8);
                                        eventCells_quatary.timeShift.ID(addEvent,2) = {experiment};
                                        eventCells_quatary.timeShift.ID(addEvent,3) = {position};
                                        eventCells_quatary.timeShift.ID(addEvent,4) = {cell};
                                        
                                        eventCells_quatary.noShift.ID(addEvent,1) = sprintf('Column_%d',index_8);
                                        eventCells_quatary.noShift.ID(addEvent,2) = {experiment};
                                        eventCells_quatary.noShift.ID(addEvent,3) = {position};
                                        eventCells_quatary.noShift.ID(addEvent,4) = {cell};
                                        index_8 = index_8 + 1;
                                        % look for next valid event or
                                        % apoptosis --> mask remaining trax
                                        break
                                    else
                                        presentEventOrder = presentEventOrder + 1;
                                    end
                            end
                        end % iCellEvents
                    end
                    
                    if ip.collect.mitosisCells_evenCells_timeShifted == 1 % mitosis at '0'
                        cellMitosis = fieldnames(cellDataStructure.(experiment).(position).(cell).mitosis);
                        for iMitosis = 1:numel(cellMitosis)
                            mitosisStatus = cellDataStructure.(experiment).(position).(cell).mitosis.(cellMitosis{iMitosis}).status;
                            switch mitosisStatus
                                case 'detected'
                                    for iTraces = 1:numel(traces)
                                        for iChannels = 1:numel(channels)
                                            mitosis_eventCells.(cellMitosis{iMitosis}).data.(traces{iTraces}).(channels{iChannels})(1:1000,index_9) = NaN;
                                            data = cellDataStructure.(experiment).(position).(cell).tsTraces.mitosis.(cellMitosis{iMitosis}).(traces{iTraces}).(channels{iChannels});
                                            mitosis_eventCells.(cellMitosis{iMitosis}).data.(traces{iTraces}).(channels{iChannels})(1:length(data),index_9) = data;
                                        end
                                    end
                                    mitosis_eventCells.(cellMitosis{iMitosis}).ID(addEvent,1) = sprintf('Column_%d',index_9);
                                    mitosis_eventCells.(cellMitosis{iMitosis}).ID(addEvent,2) = {experiment};
                                    mitosis_eventCells.(cellMitosis{iMitosis}).ID(addEvent,3) = {position};
                                    mitosis_eventCells.(cellMitosis{iMitosis}).ID(addEvent,4) = {cell};
                                    index_9 = index_9 + 1;
                            end
                            
                        end
                    end
                    
                    if ip.collect.nCellsWithEvents == 1
                        cellMitosis = fieldnames(cellDataStructure.(experiment).(position).(cell).mitosis);
                        for iMitosis = 1:numel(cellMitosis)
                            mitosisStatus = cellDataStructure.(experiment).(position).(cell).mitosis.(cellMitosis{iMitosis}).status;
                            switch mitosisStatus
                                case 'detected'
                                    %%%%%%%%%%%%%%????????
                                    mitosisFrame = cellDataStructure.(experiment).(position).(cell).mitosis.(cellMitosis{iMitosis}).frame;
                            end
                        end
                    end
                    
                    if ip.collect.apoptoticEventCells == 1 %apoptosis at '0'
                        apoptosisStatus = fieldnames(cellDataStructure.(experiment).(position).(cell).apoptosis.status);
                        switch apoptosisStatus
                            case 'detected'
                                for iTraces = 1:numel(traces)
                                    for iChannels = 1:numel(channels)
                                        apoptosis_eventCells.data.(traces{iTraces}).(channels{iChannels})(1:1000,index_12) = NaN;
                                        data = cellDataStructure.(experiment).(position).(cell).tsTraces.apoptosis.(traces{iTraces}).(channels{iChannels});
                                        apoptosis_eventCells.data.(traces{iTraces}).(channels{iChannels})(1:length(data),index_12) = data;
                                    end
                                end
                                apoptosis_eventCells.ID(addEvent,1) = sprintf('Column_%d',index_12);
                                apoptosis_eventCells.ID(addEvent,2) = {experiment};
                                apoptosis_eventCells.ID(addEvent,3) = {position};
                                apoptosis_eventCells.ID(addEvent,4) = {cell};
                                index_12 = index_12 + 1;
                        end
                    end
            end
            
        end
        
        cd ../..
    end
    clear tsDirs
    cd ../
end
fprintf('\nDone.\n')
end

