function [cellDataStructure] = detectMitosisApoptosis(ip,cellDataStructure)
%% Detects apoptosis or mitosis from nucleus size data

mitosisOrder = {'first' 'second' 'third' 'fourth' 'fifth' 'sixth' 'seventh' 'eighth' 'nineth' 'tenth'};
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
    
    data = cellDataStructure.(experiment).(position).(cellID).mfNucleusArea;
    nFrames = length(data);
    iCellMitosis = 1;  
    apoptosisDetected = 0;
    iFrame = 1;
    while iFrame <= nFrames
        if data(iFrame) < 1100
            
            % Detect apoptosis            
            if iFrame + 25 > nFrames
                apoptosisDetected = 1;
            elseif isnan(data(iFrame+25)) == 1
                apoptosisDetected = 1;
            elseif mean(data(iFrame+20:iFrame+25)) < 1300
                apoptosisDetected = 1;
            end
                   
            if apoptosisDetected == 1   
                cellDataStructure.(experiment).(position).(cellID).apoptosis.status = 'detected';
                cellDataStructure.(experiment).(position).(cellID).apoptosis.frame = iFrame;
                
                if isfield(cellDataStructure.(experiment).(position).(cellID),'mitosis') == 0
                   cellDataStructure.(experiment).(position).(cellID).mitosis.(mitosisOrder{iCellMitosis}).status = 'not detected';
                   cellDataStructure.(experiment).(position).(cellID).mitosis.(mitosisOrder{iCellMitosis}).frame = NaN;
                end
                
                break
                    
            % If not apoptosis, mitosis detected                
            else
                cellDataStructure.(experiment).(position).(cellID).mitosis.(mitosisOrder{iCellMitosis}).status = 'detected';
                cellDataStructure.(experiment).(position).(cellID).mitosis.(mitosisOrder{iCellMitosis}).frame = iFrame;
                iCellMitosis = iCellMitosis + 1; 
                
                cellDataStructure.(experiment).(position).(cellID).apoptosis.status = 'not detected';
                cellDataStructure.(experiment).(position).(cellID).apoptosis.frame = NaN;
                
                iFrame = iFrame + 25;
                continue
            end        
        end
        iFrame = iFrame + 1;
    end
    
    if isfield(cellDataStructure.(experiment).(position).(cellID),'mitosis') == 0
        cellDataStructure.(experiment).(position).(cellID).mitosis.(mitosisOrder{iCellMitosis}).status = 'not detected';
        cellDataStructure.(experiment).(position).(cellID).mitosis.(mitosisOrder{iCellMitosis}).frame = NaN;
    end
    
    if isfield(cellDataStructure.(experiment).(position).(cellID),'apoptosis') == 0
        cellDataStructure.(experiment).(position).(cellID).apoptosis.status = 'not detected';
        cellDataStructure.(experiment).(position).(cellID).apoptosis.frame = NaN;
    end
end
cd (ip.fdp)
end