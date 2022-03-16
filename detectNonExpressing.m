function [cellDataStructure] = detectNonExpressing(ip,cellDataStructure)
%% Detect cells not expressing eGFP marker
expressionCutOff = 0.001;
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
    
    data = cellDataStructure.(experiment).(position).(cellID).mfBcTraces.eGFP;
    firstFrame = cellDataStructure.(experiment).(position).(cellID).cellTrackingStarted;
    
    if isnan(firstFrame) == 1
        cellDataStructure.(experiment).(position).(cellID).eGFPExpressionLevel = NaN; 
    elseif data(firstFrame) < expressionCutOff
        cellDataStructure.(experiment).(position).(cellID).eGFPExpressionLevel = 'subthreshold';
    else
        cellDataStructure.(experiment).(position).(cellID).eGFPExpressionLevel = 'adequate';
    end

end
cd (ip.fdp)
end