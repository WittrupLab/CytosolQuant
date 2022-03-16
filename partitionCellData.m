function [cellDataStructure] = partitionCellData(ip)
%% Rearranging CellProfiler tracking data
experimentFolder = [ip.dataDir filesep ip.listAcquisitions{ip.iAcquisitions,1}];

columIndex = [12,73; 13,72; 13,73; 26,1; 21,1; 22,1];
sheetFolderNames = {...
    'identifyNuc' 'cellMaskedsiRNA' 'backgroundIntensity' ;
    'identifyNuc' 'identifyNuc' 'backgroundIntensity';
    'identifyNuc' 'cellFullsiRNA' 'backgroundIntensity';
    'identifyNuc' 'identifyNuc' 'backgroundIntensity';
    'identifyNuc' 'identifyNuc' 'backgroundIntensity';
    'identifyNuc' 'identifyNuc' 'backgroundIntensity'};

for ll = 1:6
    
    folder = [experimentFolder filesep ip.listAcquisitions{ip.iAcquisitions,2} filesep 'cellprofilerOutput' filesep 'dataSheet' filesep char(sheetFolderNames{ll,1})];
    cd (folder)
    clear dir
    dirData = dir; % Get the data for the current directory
    dirData(strncmp({dirData.name}, '.', 1)) = []; % Remove false files
    dirData(strncmp({dirData.name}, '~', 1)) = [];
    fileName = char({dirData.name}');
    sheet = table2array(readtable(fileName));
    sheet(1:end, 1) = sheet(1:end, 1)-sheet(1,1)+1;
    frames = sheet(1:end, 1);
    nFrames = max(frames);
    trackID = sheet(1:end, 65);
    
    folder = [experimentFolder filesep ip.listAcquisitions{ip.iAcquisitions,2} filesep 'cellprofilerOutput' filesep 'dataSheet' filesep char(sheetFolderNames{ll,2})];
    cd (folder)
    clear dir
    dirData = dir; % Get the data for the current directory
    dirData(strncmp({dirData.name}, '.', 1)) = []; % Remove false files
    dirData(strncmp({dirData.name}, '~', 1)) = [];
    fileName = char({dirData.name}');
    sheet2 = table2array(readtable(fileName));
    medianintensity = sheet2(1:end, columIndex(ll,1));
    
    folder = [experimentFolder filesep ip.listAcquisitions{ip.iAcquisitions,2} filesep 'cellprofilerOutput' filesep 'dataSheet' filesep char(sheetFolderNames{ll,3})];
    cd (folder)
    clear dir
    dirData = dir; % Get the data for the current directory
    dirData(strncmp({dirData.name}, '.', 1)) = []; % Remove false files
    dirData(strncmp({dirData.name}, '~', 1)) = [];
    fileName = char({dirData.name}');
    sheet3 = readtable(fileName);
    backgroundintensity = table2array(sheet3(1:end, columIndex(ll,2)));
    
  resampledValues = NaN(nFrames, max(trackID));
    
    for ii = 1:size(sheet,1)
        
        if isnan(trackID(ii)) == 1
            continue
        end
        
        tempID_LAP = trackID(ii);
        tempID_frame = sheet(ii, 2);
        tempframe = sheet(ii, 1);
        
        if isnan(tempID_LAP) == 1
            continue
        end
        
        for jj = 1:size(sheet2,1)
            
            if tempID_frame == sheet(jj, 2) && ...
                    tempframe == sheet(jj, 1)
                
                resampledValues(tempframe, tempID_LAP) = medianintensity(jj);
            end
        end
    end
    
    nanRef = resampledValues;
    nanRef(~isnan(nanRef)) = 1;
    if ll == 1
        siRNA = resampledValues; 
        BCsiRNA = resampledValues-backgroundintensity;
        MFBCsiRNA = movmedian(BCsiRNA,5,'omitnan') .* nanRef;
        MFsiRNA = movmedian(siRNA,5,'omitnan') .* nanRef;
        bg_siRNA = backgroundintensity;
    elseif ll == 2
        d1eGFP = resampledValues; 
        BCd1eGFP = resampledValues-backgroundintensity;
        MFBCd1eGFP = movmedian(BCd1eGFP,5,'omitnan') .* nanRef;
        MFd1eGFP = movmedian(d1eGFP,5,'omitnan') .* nanRef;
        bg_d1eGFP = backgroundintensity;
    elseif ll == 3
        siRNAF = resampledValues;
        BCsiRNAF = resampledValues-backgroundintensity;
        MFBCsiRNAF = movmedian(BCsiRNAF,5,'omitnan') .* nanRef;
        MFsiRNAF = movmedian(siRNAF,5,'omitnan') .* nanRef;
        bg_siRNAF = backgroundintensity;
    elseif ll == 4
        nucleusArea = resampledValues;
        MFnucleusArea = movmedian(nucleusArea,5,'omitnan') .* nanRef;
    elseif ll == 5
        xCellCoordinates = resampledValues;
    elseif ll == 6
        yCellCoordinates = resampledValues;
    end
    
end

folder = [experimentFolder filesep ip.listAcquisitions{ip.iAcquisitions,2} filesep 'matlabOutput'];
cd (folder)
filename = char(ip.listAcquisitions{ip.iAcquisitions,2});
if (filename(length(filename)) == '#') == 1
    filename = filename(1:(length(filename)-1));
end

load(sprintf('%s%s', filename, '_defineControls'));

for iCells = 1:max(trackID)
    experiment = ip.listAcquisitions{ip.iAcquisitions,1};
    position = ip.listAcquisitions{ip.iAcquisitions,2};
    if (position(length(position)) == '#') == 1
        position = position(1:(length(position)-1));
    end
    POS = position(end-3:end);
    controlPositionStatus = controlPositionLogicTable.(POS)(1);
    
    if iCells < 10
        cellID = sprintf('cell_00%d',iCells);    
    elseif iCells < 100
        cellID = sprintf('cell_0%d',iCells);   
    else
        cellID = sprintf('cell_%d',iCells);
    end
    
    cellDataStructure.(experiment).(position).(cellID).rawTraces.eGFP = d1eGFP(:,iCells);
    cellDataStructure.(experiment).(position).(cellID).rawTraces.siRNA = siRNA(:,iCells);
    cellDataStructure.(experiment).(position).(cellID).rawTraces.siRNAF = siRNAF(:,iCells);
    
    cellDataStructure.(experiment).(position).(cellID).bcTraces.eGFP = BCd1eGFP(:,iCells);
    cellDataStructure.(experiment).(position).(cellID).bcTraces.siRNA = BCsiRNA(:,iCells);
    cellDataStructure.(experiment).(position).(cellID).bcTraces.siRNAF = BCsiRNAF(:,iCells);
    
    cellDataStructure.(experiment).(position).(cellID).mfTraces.eGFP = MFd1eGFP(:,iCells);
    cellDataStructure.(experiment).(position).(cellID).mfTraces.siRNA = MFsiRNA(:,iCells);
    cellDataStructure.(experiment).(position).(cellID).mfTraces.siRNAF = MFsiRNAF(:,iCells);
    
    cellDataStructure.(experiment).(position).(cellID).mfBcTraces.eGFP = MFBCd1eGFP(:,iCells);
    cellDataStructure.(experiment).(position).(cellID).mfBcTraces.siRNA = MFBCsiRNA(:,iCells);
    cellDataStructure.(experiment).(position).(cellID).mfBcTraces.siRNAF = MFBCsiRNAF(:,iCells);
    
    cellDataStructure.(experiment).(position).(cellID).rawNucleusArea = nucleusArea(:,iCells);
    cellDataStructure.(experiment).(position).(cellID).mfNucleusArea =  MFnucleusArea(:,iCells);
    cellDataStructure.(experiment).(position).(cellID).backgroundIntensity.eGFP = bg_d1eGFP;
    cellDataStructure.(experiment).(position).(cellID).backgroundIntensity.siRNA = bg_siRNA;
    cellDataStructure.(experiment).(position).(cellID).backgroundIntensity.siRNAF = bg_siRNAF;
    
    if controlPositionStatus == 0
        cellDataStructure.(experiment).(position).(cellID).controlPosition = 'no';
    elseif controlPositionStatus == 1
        cellDataStructure.(experiment).(position).(cellID).controlPosition = 'yes';
    end
        
    
    P = {'xCoordinate' 'yCoordinate'};
    coordinateArray(:,1) = xCellCoordinates(:,iCells);
    coordinateArray(:,2) = yCellCoordinates(:,iCells);
    coordinateTable = array2table(coordinateArray);
    coordinateTable.Properties.VariableNames = P;
    cellDataStructure.(experiment).(position).(cellID).trackingCoordinates = coordinateTable;
end
cd (ip.fdp)
end

