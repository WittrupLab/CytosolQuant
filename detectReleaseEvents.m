function [cellDataStructure] = detectReleaseEvents(ip,cellDataStructure)
%% Detects siRNA release events 

%Outputs array of cellIDs and frame of first event defined as
%increase in signal intensity above 3 stdev from input matrix organized as
%frames in 1st dim and cellIDs in 2nd dim.
experiment = ip.listAcquisitions{ip.iAcquisitions,1};
position = ip.listAcquisitions{ip.iAcquisitions,2};
if (position(length(position)) == '#') == 1
   position = position(1:(length(position)-1));
end
nCells = length(fieldnames(cellDataStructure.(experiment).(position)));

for iCells = 1:nCells
    proceed = 0;
    eventDetected = 0;
    
    if iCells < 10
        cellID = sprintf('cell_00%d',iCells);    
    elseif iCells < 100
        cellID = sprintf('cell_0%d',iCells);   
    else
        cellID = sprintf('cell_%d',iCells);
    end
    
    
    data = cellDataStructure.(experiment).(position).(cellID).rawTraces.eGFP;
    nFrames = length(data);
    
    % Find first tracked frame in raw data traces (not NaN)
    for iFrame = 1:nFrames
        if isnan(data(iFrame)) == 0
           cellDataStructure.(experiment).(position).(cellID).cellTrackingStarted = iFrame;
           
           if iFrame > 1
               cellDataStructure.(experiment).(position).(cellID).mfTraces.eGFP(1:iFrame-1) = NaN;
               cellDataStructure.(experiment).(position).(cellID).mfTraces.siRNA(1:iFrame-1) = NaN;
               cellDataStructure.(experiment).(position).(cellID).mfTraces.siRNAF(1:iFrame-1) = NaN;
               cellDataStructure.(experiment).(position).(cellID).mfBcTraces.eGFP(1:iFrame-1) = NaN;
               cellDataStructure.(experiment).(position).(cellID).mfBcTraces.siRNA(1:iFrame-1) = NaN;
               cellDataStructure.(experiment).(position).(cellID).mfBcTraces.siRNAF(1:iFrame-1) = NaN;
           end   
           break
           
        elseif iFrame == nFrames
            cellDataStructure.(experiment).(position).(cellID).cellTrackingStarted = NaN; 
        end
    end    
    
    % Find last tracked frame in raw data traces (not NaN)
    for iFrame = 1:nFrames
        revFrame = nFrames - iFrame + 1;
        if isnan(data(revFrame)) == 0
           cellDataStructure.(experiment).(position).(cellID).cellTrackingEnded = revFrame;
           
           if revFrame < nFrames
               cellDataStructure.(experiment).(position).(cellID).mfTraces.eGFP(revFrame+1:nFrames) = NaN;
               cellDataStructure.(experiment).(position).(cellID).mfTraces.siRNA(revFrame+1:nFrames) = NaN;
               cellDataStructure.(experiment).(position).(cellID).mfTraces.siRNAF(revFrame+1:nFrames) = NaN;
               cellDataStructure.(experiment).(position).(cellID).mfBcTraces.eGFP(revFrame+1:nFrames) = NaN;
               cellDataStructure.(experiment).(position).(cellID).mfBcTraces.siRNA(revFrame+1:nFrames) = NaN;
               cellDataStructure.(experiment).(position).(cellID).mfBcTraces.siRNAF(revFrame+1:nFrames) = NaN;
           end
           
           break
           
        elseif iFrame == nFrames
            cellDataStructure.(experiment).(position).(cellID).cellTrackingEnded = NaN;
        end
    end
    
    % Calculate duration of track 
    firstFrame = cellDataStructure.(experiment).(position).(cellID).cellTrackingStarted;
    endFrame = cellDataStructure.(experiment).(position).(cellID).cellTrackingEnded;
    trackDuration = endFrame - firstFrame + 1;
    cellDataStructure.(experiment).(position).(cellID).trackDuration = trackDuration;
    
    % Find gaps (NaNs) in cell track
    if isnan(firstFrame) == 0 && isnan(endFrame) == 0
        trackingGaps = sum(isnan(data(firstFrame:endFrame)));
        cellDataStructure.(experiment).(position).(cellID).trackingGaps = trackingGaps;
    else
        cellDataStructure.(experiment).(position).(cellID).trackingGaps = NaN;
    end
    
    data = cellDataStructure.(experiment).(position).(cellID).mfBcTraces.siRNAF;
    
    % Find first NonNaN value in median filtered data traces
    % Proceed to detect release events
    firstNonNaN = 1;
    while isnan(data(firstNonNaN)) == 1
        if firstNonNaN < nFrames - 10
        firstNonNaN = firstNonNaN + 1;
        else
            cellDataStructure.(experiment).(position).(cellID).allEvents.first.status = 'not detected';
            cellDataStructure.(experiment).(position).(cellID).allEvents.first.frame = NaN;
            cellDataStructure.(experiment).(position).(cellID).allEvents.first.validationStatus = NaN;
            proceed = 1;
            break
        end
    end
        
    if proceed == 1
        continue
    end
    
    eventOrder = {'first' 'second' 'third' 'fourth' 'fifth' 'sixth' 'seventh' 'eighth' 'nineth' 'tenth'};
    iCellEvents = 1;
    if firstNonNaN == 1
       evalFrame = firstNonNaN+3;
    else
       evalFrame = firstNonNaN+5;
    end
    
    while evalFrame < nFrames - 10
        if firstNonNaN > evalFrame - 10
           stdev = std(data(firstNonNaN:evalFrame-1));
        else
           stdev = std(data(evalFrame-10:evalFrame-1));
        end
        
          preReleaseMean = mean(data(evalFrame-3:evalFrame-1));
         
          if all(data(evalFrame:evalFrame+5)>=max(stdev*3+preReleaseMean,0.00003+preReleaseMean))
         
             cellDataStructure.(experiment).(position).(cellID).allEvents.(eventOrder{iCellEvents}).status = 'detected';
             cellDataStructure.(experiment).(position).(cellID).allEvents.(eventOrder{iCellEvents}).frame = evalFrame;
             cellDataStructure.(experiment).(position).(cellID).allEvents.(eventOrder{iCellEvents}).validationStatus = 'not validated';
                 
             evalFrame = evalFrame+9;
             iCellEvents = iCellEvents + 1;
             eventDetected = 1;
          else
             evalFrame = evalFrame+1;
          end
    end 
        
    if eventDetected == 0
       cellDataStructure.(experiment).(position).(cellID).allEvents.(eventOrder{iCellEvents}).status = 'not detected';
       cellDataStructure.(experiment).(position).(cellID).allEvents.(eventOrder{iCellEvents}).frame = NaN;
       cellDataStructure.(experiment).(position).(cellID).allEvents.(eventOrder{iCellEvents}).validationStatus = NaN;
    end
    cellDataStructure.(experiment).(position).(cellID).disqualifyCell = 'no';
            
end
cd (ip.fdp)
end