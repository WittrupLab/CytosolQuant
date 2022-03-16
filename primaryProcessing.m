function [] = primaryProcessing(ip)
%% Primary processing (batch-loop, only performed once per experiment is fixed parameters)
fprintf('\nPrimary processing initiated...\n\n');
cd (ip.fdp)
ip.listAcquisitions = listAcquisitions(ip);

for iAcquisitions = 1:size(ip.listAcquisitions,1)
    ip.iAcquisitions = iAcquisitions;
    fprintf('Processing %s\n',ip.listAcquisitions{iAcquisitions,2});
    position = ip.listAcquisitions{ip.iAcquisitions,2};
    if (position(length(position)) == '#') == 1
        position = position(1:(length(position)-1));
    end
    fprintf('%s\n', position);
    
    
    % 1.Sorting of raw exported images and data sheets in folder organization
    cd (ip.fdp)
    copyImages(ip)
    
    cd (ip.fdp)
    sortResampleImages(ip)
    
    % 2.Partition of collective datasheets in single cell structured array
    cd (ip.fdp)
    cellDataStructure = partitionCellData(ip);
    
    % 3.Detection of release events defining event cells and non-event cells
    % Determine duration of track, first and last tracked frame, track gaps
    cellDataStructure = detectReleaseEvents(ip,cellDataStructure);
    
    % 4.Detect mitosis and apoptosis
    cellDataStructure = detectMitosisApoptosis(ip,cellDataStructure);
    
    % 5.Detect non-expressing cells
    cellDataStructure = detectNonExpressing(ip,cellDataStructure);
    
    % 8.Perform time-shift of release event and mitosis
    cd (ip.fdp)
    cellDataStructure = timeShiftTraces(ip,cellDataStructure);
    
    % Save cell data structure in position folder
    folder = [ip.dataDir filesep ip.listAcquisitions{ip.iAcquisitions,1} filesep ip.listAcquisitions{ip.iAcquisitions,2} filesep 'matlabOutput'];
    cd (folder)
    filename = char(sprintf('%s%s', position,'_cellDataStructure'));
    save(char(filename),'cellDataStructure')
    
    cd (ip.fdp)
end
fprintf('\nPrimary processing completed.\n');
end