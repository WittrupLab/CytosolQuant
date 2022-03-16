function [cellDataStructure] = timeShiftTraces(ip,cellDataStructure)
%% Shifts data in T-dimension so that the time-point of all events (or mitosis or apoptosis) align

eventOrder = {'first' 'second' 'third' 'fourth' 'fifth' 'sixth' 'seventh' 'eighth' 'nineth' 'tenth'};
mitosisOrder = eventOrder;
experiment = ip.listAcquisitions{ip.iAcquisitions,1};
position = ip.listAcquisitions{ip.iAcquisitions,2};
if (position(length(position)) == '#') == 1
    position = position(1:(length(position)-1));
end
nCells = length(fieldnames(cellDataStructure.(experiment).(position)));

for iCells = 1:nCells
    
    if iCells < 10
        cellID = sprintf('cell_00%d',iCells);
    elseif iCells < 100
        cellID = sprintf('cell_0%d',iCells);
    else
        cellID = sprintf('cell_%d',iCells);
    end
    
    firstFrame = cellDataStructure.(experiment).(position).(cellID).cellTrackingStarted;
    endFrame = cellDataStructure.(experiment).(position).(cellID).cellTrackingEnded;
    allCellEvents = fieldnames(cellDataStructure.(experiment).(position).(cellID).allEvents);
    allDataTypes = {'rawTraces' 'mfTraces' 'mfBcTraces' 'bcTraces'};
    allParameters = {'eGFP' 'siRNA' 'siRNAF'};
    
    % Time shift according to event frame
    for iCellEvents = 1:length(allCellEvents)
        eventStatus = cellDataStructure.(experiment).(position).(cellID).allEvents.(eventOrder{iCellEvents}).status;
        eventFrame = cellDataStructure.(experiment).(position).(cellID).allEvents.(eventOrder{iCellEvents}).frame;
        for iDataTypes = 1:length(allDataTypes)
            for iParameters = 1:length(allParameters)
                switch eventStatus
                    case 'detected'
                        data = cellDataStructure.(experiment).(position).(cellID).(allDataTypes{iDataTypes}).(allParameters{iParameters});
                        shiftedData = NaN(1000,1);
                        shiftedData(500-(eventFrame-firstFrame):500+(endFrame-eventFrame)) = data(firstFrame:endFrame);
                        cellDataStructure.(experiment).(position).(cellID).tsTraces.events.(eventOrder{iCellEvents}).(allDataTypes{iDataTypes}).(allParameters{iParameters}) = shiftedData;
                    case {'not detected', 'ab inito'}
                        shiftedData = NaN(1000,1);
                        cellDataStructure.(experiment).(position).(cellID).tsTraces.events.(eventOrder{iCellEvents}).(allDataTypes{iDataTypes}).(allParameters{iParameters}) = shiftedData;
                end
            end
        end
    end
    
    % Time shift according to mitosis detection
    allMitosis = fieldnames(cellDataStructure.(experiment).(position).(cellID).mitosis);
    for iMitosis = 1:length(allMitosis)
        mitosisStatus = cellDataStructure.(experiment).(position).(cellID).mitosis.(mitosisOrder{iMitosis}).status;
        mitosisFrame = cellDataStructure.(experiment).(position).(cellID).mitosis.(mitosisOrder{iMitosis}).frame;
        for iDataTypes = 1:length(allDataTypes)
            for iParameters = 1:length(allParameters)
                switch mitosisStatus
                    case 'detected'
                        data = cellDataStructure.(experiment).(position).(cellID).(allDataTypes{iDataTypes}).(allParameters{iParameters});
                        shiftedData = NaN(1000,1);
                        shiftedData(500-(mitosisFrame-firstFrame):500+(endFrame-mitosisFrame)) = data(firstFrame:endFrame);
                        cellDataStructure.(experiment).(position).(cellID).tsTraces.mitosis.(mitosisOrder{iMitosis}).(allDataTypes{iDataTypes}).(allParameters{iParameters}) = shiftedData;
                    case 'not detected'
                        shiftedData = NaN(1000,1);
                        cellDataStructure.(experiment).(position).(cellID).tsTraces.mitosis.(mitosisOrder{iMitosis}).(allDataTypes{iDataTypes}).(allParameters{iParameters}) = shiftedData;
                end
            end
        end
    end
    
    % Time shift according to apoptosis detection
    apoptosisStatus = cellDataStructure.(experiment).(position).(cellID).apoptosis.status;
    apoptosisFrame = cellDataStructure.(experiment).(position).(cellID).apoptosis.frame;
    for iDataTypes = 1:length(allDataTypes)
        for iParameters = 1:length(allParameters)
            switch apoptosisStatus
                case 'detected'
                    data = cellDataStructure.(experiment).(position).(cellID).(allDataTypes{iDataTypes}).(allParameters{iParameters});
                    shiftedData = NaN(1000,1);
                    shiftedData(500-(apoptosisFrame-firstFrame):500+(endFrame-apoptosisFrame)) = data(firstFrame:endFrame);
                    cellDataStructure.(experiment).(position).(cellID).tsTraces.apoptosis.(allDataTypes{iDataTypes}).(allParameters{iParameters}) = shiftedData;
                case 'not detected'
                    shiftedData = NaN(1000,1);
                    cellDataStructure.(experiment).(position).(cellID).tsTraces.apoptosis.(allDataTypes{iDataTypes}).(allParameters{iParameters}) = shiftedData;
            end
        end
    end
  
end
cd (ip.fdp)
end