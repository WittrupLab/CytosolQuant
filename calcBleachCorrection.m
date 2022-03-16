function [] = calcBleachCorrection( ip )
%% Calculate and correct bleaching in GFP channel
cd (ip.fdp)
ip.listAcquisitions = listAcquisitions(ip);

for iAcquisitions = 1:size(ip.listAcquisitions,1)
    ip.iAcquisitions = iAcquisitions;
    fprintf('Processing %s\n',ip.listAcquisitions{iAcquisitions,2});
    
    folder = [ip.dataDir filesep ip.listAcquisitions{ip.iAcquisitions,1} filesep ip.listAcquisitions{ip.iAcquisitions,2} filesep 'matlabOutput'];
    cd (folder)
    
    posName = char(ip.listAcquisitions{ip.iAcquisitions,2});
    if (posName(length(posName)) == '#') == 1
        posName = posName(1:(length(posName)-1));
    end
    
    filename = char(sprintf('%s%s', posName,'_intensityArray'));
    
    experiment = ip.listAcquisitions{ip.iAcquisitions,1};
    position = ip.listAcquisitions{ip.iAcquisitions,2};
    if (position(length(position)) == '#') == 1
        position = position(1:(length(position)-1));
    end
    
    data.(char(experiment)).(char(position)) = load(filename,'MFsiRNA','MFsiRNAF','MFd1eGFP','MFnucleusArea');
    
end



% Remove dying cells
for iAcquisitions = 1:size(ip.listAcquisitions,1)
    ip.iAcquisitions = iAcquisitions;
    fprintf('Processing %s\n',ip.listAcquisitions{iAcquisitions,2});
    
    experiment = ip.listAcquisitions{ip.iAcquisitions,1};
    position = ip.listAcquisitions{ip.iAcquisitions,2};
    if (position(length(position)) == '#') == 1
        position = position(1:(length(position)-1));
    end
    
    nFrames = size(data.(char(experiment)).(char(position)).MFsiRNA, 1);
    nCells = size(data.(char(experiment)).(char(position)).MFsiRNA, 2);
    data.(char(experiment)).(char(position)).mitosis = zeros(0);
    
    for i=1:nCells
        j=1;
        while j<nFrames+1
            if data.(char(experiment)).(char(position)).MFnucleusArea(j,i) < 1100
                [i,j];
                dead = false;
                if j+25 > nFrames
                    dead = true;
                else
                    if not(mean(data.(char(experiment)).(char(position)).MFnucleusArea(j+20:j+25, i)) > 1300)
                        dead = true;
                    end
                end
                
                if dead
                    firstFrameToRemove = max(1, j-20);
                    data.(char(experiment)).(char(position)).MFsiRNA(firstFrameToRemove:end, i) = NaN;
                    data.(char(experiment)).(char(position)).MFsiRNAF(firstFrameToRemove:end, i) = NaN;
                    data.(char(experiment)).(char(position)).MFd1eGFP(firstFrameToRemove:end, i) = NaN;
                    j=nFrames;
                end
                if not(dead)
                    data.(char(experiment)).(char(position)).mitosis = vertcat(data.(char(experiment)).(char(position)).mitosis, [i,j]);
                end
                j=j+25;
            end
            j=j+1;
        end
    end
end

% Remove non-expressing cells
expressionCutOff = 0.001;

for iAcquisitions = 1:size(ip.listAcquisitions,1)
    ip.iAcquisitions = iAcquisitions;
    fprintf('Processing %s\n',ip.listAcquisitions{iAcquisitions,2});
    
    experiment = ip.listAcquisitions{ip.iAcquisitions,1};
    position = ip.listAcquisitions{ip.iAcquisitions,2};
    if (position(length(position)) == '#') == 1
        position = position(1:(length(position)-1));
    end
    
    nCells = size(data.(char(experiment)).(char(position)).MFsiRNA, 2);
    cd (ip.fdp)
    for i=1:nCells
        if firstNonNaN(data.(char(experiment)).(char(position)).MFd1eGFP(1:end, i)) < expressionCutOff
            data.(char(experiment)).(char(position)).MFsiRNA(1:end, i) = NaN;
            data.(char(experiment)).(char(position)).MFsiRNAF(1:end, i) = NaN;
            data.(char(experiment)).(char(position)).MFd1eGFP(1:end, i) = NaN;
        end
    end
    
end


% Read data from multiple positions and shift data in time to first frame
% preallocate bleachCorrData structured array
for iAcquisitions = 1:size(ip.listAcquisitions,1)
    ip.iAcquisitions = iAcquisitions;
    fprintf('Processing %s\n',ip.listAcquisitions{iAcquisitions,2});
    
    experiment = ip.listAcquisitions{ip.iAcquisitions,1};
    bleachCorrData.(char(experiment)) = [];
end
for iAcquisitions = 1:size(ip.listAcquisitions,1)
    ip.iAcquisitions = iAcquisitions;
    fprintf('Processing %s\n',ip.listAcquisitions{iAcquisitions,2});
    
    experiment = ip.listAcquisitions{ip.iAcquisitions,1};
    position = ip.listAcquisitions{ip.iAcquisitions,2};
    if (position(length(position)) == '#') == 1
        position = position(1:(length(position)-1));
    end
    
    if iAcquisitions == 1
        bleachCorrData.(char(experiment)) = [];
    end
    
    nFrames = size(data.(char(experiment)).(char(position)).MFsiRNA, 1);
    
    %  after event at nFrames+1 % DETECT RELEASE EVENTS
    cd (ip.fdp)
    data.(char(experiment)).(char(position)).allEventFrames = detectReleaseEvents(data.(char(experiment)).(char(position)).MFsiRNAF);
    cd (ip.fdp)
    data.(char(experiment)).(char(position)).allEventFrames = excludeFalseEvents(ip, data.(char(experiment)).(char(position)).allEventFrames);
    
    % Create eventCellIndex key
    nEvents = 0;
    eventCells = data.(char(experiment)).(char(position)).allEventFrames(1, 1:end) > 0;
    for hh = 1:length(eventCells)
        if eventCells(hh) == 1
            data.(char(experiment)).(char(position)).eventCellIndexKey(nEvents+1) = hh;
            nEvents = nEvents + 1;
        end
    end
    
    % Remove cells that are not traced from start of acquisition
    nFrames_treshold = 3;
    nCells = size(data.(char(experiment)).(char(position)).MFsiRNA, 2);
    for i=1:nCells
        if sum(isnan(data.(char(experiment)).(char(position)).MFsiRNA(1:nFrames_treshold,i))) == 10
            data.(char(experiment)).(char(position)).MFsiRNA(1:end, i) = NaN;
            data.(char(experiment)).(char(position)).MFsiRNAF(1:end, i) = NaN;
            data.(char(experiment)).(char(position)).MFd1eGFP(1:end, i) = NaN;
            data.(char(experiment)).(char(position)).MFnucleusArea(1:end, i) = NaN;
            % remove cells from mitosis cell list
            isMitosis = i ~= data.(char(experiment)).(char(position)).mitosis(:,1);
            data.(char(experiment)).(char(position)).mitosis = data.(char(experiment)).(char(position)).mitosis(isMitosis,:);
        end
    end
    
    eventCells = data.(char(experiment)).(char(position)).allEventFrames(1, 1:end) > 0;
    nonEventCells = isnan(data.(char(experiment)).(char(position)).allEventFrames(1, 1:end));
    priorEventCells = data.(char(experiment)).(char(position)).allEventFrames(1, 1:end) < 0;
    secondaryEventFrames = data.(char(experiment)).(char(position)).allEventFrames(2, eventCells);
    
    %define last frame to be analysed for cells having a single event, i.e.
    %last frame or secondary event frame
    nEventCells = sum(eventCells);
    
    lastSingleEventFrame = min(nFrames, secondaryEventFrames);
    
    data.(char(experiment)).(char(position)).siRNAEventCells = data.(char(experiment)).(char(position)).MFsiRNA(1:end, eventCells);
    for i=1:nEventCells
        data.(char(experiment)).(char(position)).siRNAEventCells(lastSingleEventFrame(i):end, i) = NaN;
    end
    data.(char(experiment)).(char(position)).eGFPEventCells = data.(char(experiment)).(char(position)).MFd1eGFP(1:end, eventCells);
    for i=1:nEventCells
        data.(char(experiment)).(char(position)).eGFPEventCells(lastSingleEventFrame(i):end, i) = NaN;
    end
    data.(char(experiment)).(char(position)).eGFPNonEventCells = data.(char(experiment)).(char(position)).MFd1eGFP(1:end, nonEventCells);
    addData = data.(char(experiment)).(char(position)).MFd1eGFP(1:end, nonEventCells);
    corrData = bleachCorrData.(char(experiment));
    bleachCorrData.(char(experiment)).MFd1eGFP(1:size(addData,1),size(corrData,2)+1:size(corrData,2)+size(addData,2)) = addData;
    
end

% Normalize for bleaching, trim away non-expressing non-event-cells
for iAcquisitions = 1:size(ip.listAcquisitions,1)
    ip.iAcquisitions = iAcquisitions;
    fprintf('Processing %s\n',ip.listAcquisitions{iAcquisitions,2});
    
    experiment = ip.listAcquisitions{ip.iAcquisitions,1};
    position = ip.listAcquisitions{ip.iAcquisitions,2};
    if (position(length(position)) == '#') == 1
        position = position(1:(length(position)-1));
    end
    
    nFrames = size(data.(char(experiment)).(char(position)).MFsiRNA, 1);
    
    cd (ip.fdp)
    eventCells = data.(char(experiment)).(char(position)).allEventFrames(1, 1:end) > 0;
    nonEventCells = isnan(data.(char(experiment)).(char(position)).allEventFrames(1, 1:end));
    
    bleachCorrData.(char(experiment)).trimedNonEventeGFP = bleachCorrData.(char(experiment)).MFd1eGFP;
    
    removeNonexpressing = 0;
    nNonEventCells = size(bleachCorrData.(char(experiment)).trimedNonEventeGFP, 2);
    for i=1:nNonEventCells
        if (nanmean(bleachCorrData.(char(experiment)).trimedNonEventeGFP(1:nFrames, i)) < 0.001)
            bleachCorrData.(char(experiment)).trimedNonEventeGFP(1:end, i) = NaN;
            removeNonexpressing = removeNonexpressing +1;
        end
    end
    bleachCorrData.(char(experiment)).removeNonexpressing = removeNonexpressing;
end

% Trim away short non-event tracks
for iAcquisitions = 1:size(ip.listAcquisitions,1)
    ip.iAcquisitions = iAcquisitions;
    fprintf('Processing %s\n',ip.listAcquisitions{iAcquisitions,2});
    
    experiment = ip.listAcquisitions{ip.iAcquisitions,1};
    position = ip.listAcquisitions{ip.iAcquisitions,2};
    if (position(length(position)) == '#') == 1
        position = position(1:(length(position)-1));
    end
    
    trackLengthNonEvent = sum(not(isnan(bleachCorrData.(char(experiment)).trimedNonEventeGFP)), 1);
    
    nNonEventCells = size(bleachCorrData.(char(experiment)).trimedNonEventeGFP, 2) - ...
        bleachCorrData.(char(experiment)).removeNonexpressing;
    for i=1:nNonEventCells
        if trackLengthNonEvent(i) < 50
            bleachCorrData.(char(experiment)).trimedNonEventeGFP(1:end, i) = NaN;   %%Remove tracks shorter than 50
        end
    end
    
    bleachCorrFactor = nanmean(bleachCorrData.(char(experiment)).trimedNonEventeGFP, 2);
    bleachCorrData.(char(experiment)).bleachCorrFactor = bleachCorrFactor ./ bleachCorrFactor(1);
end

% Save
for iAcquisitions = 1:size(ip.listAcquisitions,1)
    ip.iAcquisitions = iAcquisitions;
    fprintf('Processing %s\n',ip.listAcquisitions{iAcquisitions,2});
    
    experiment = ip.listAcquisitions{ip.iAcquisitions,1};
    position = ip.listAcquisitions{ip.iAcquisitions,2};
    if (position(length(position)) == '#') == 1
        position = position(1:(length(position)-1));
    end
    
    detectedReleaseEvents = data.(char(experiment)).(char(position));
    folder = [ip.dataDir filesep ip.listAcquisitions{ip.iAcquisitions,1} filesep ip.listAcquisitions{ip.iAcquisitions,2} filesep 'matlabOutput'];
    cd (folder)
    
    filename = char(sprintf('%s%s', position,'_bleachCorrectionData'));
    save(char(filename),'bleachCorrData')
end


