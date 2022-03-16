function [] = addValidationData( ip, readpath )
%% Partition manual tracks
fprintf('\nAdding validation status...\n\n')
cd (readpath)
load('panelMetadata.mat');

allChannels = ip.allChannels.(char(ip.allExperiments(1)));
ip.nChannels = length(allChannels);

if ip.noFalseEvents == 0
    
    % Read the manual Tracks file from ImageJ
    table=readtable('falseEvents.csv');
    table.Type = [];
    table.Color = [];
    table.Fill = [];
    
    % slice no is dependent on which channel is selected in imagej while
    % tracking (as well as the number of z-planes and channels). The routine is
    % 1) perform object tracking  choossing the first channel 2) and prior to
    % automatic object tracking
    % Thus, there can be 2 or 3 de facto channels.
    
    % table(:,1) listindex
    % table(:,2) track no
    % table(:,3) slice no
    % table(:,4) X coordinate
    % table(:,5) Y coordinate
    % table(:,6) distance
    % table(:,7) velocity
    % table(:,8) pixel value
    
    fprintf('\n\nFalse or disqualified events\n')
    for ii = 1:size(table,1)
        % convert slice no to relative time variable
        fprintf('.');
        
        xVar = table2array(table(ii,4));
        yVar = table2array(table(ii,5));
        
        for jj = 1:size(panelMetadata,1)
            
            xLim = panelMetadata{jj,10};
            yLim = panelMetadata{jj,11};
            xDim = panelMetadata{jj,5};
            yDim = panelMetadata{jj,6};
            
            % evaluate if the coordinates belong to the current ROI
            if xVar >= xLim && xVar <= xLim+xDim   && ...
                    yVar >= yLim && yVar <= yLim+yDim
                
                exp = panelMetadata{jj,1};
                ts = panelMetadata{jj,2};
                cell = panelMetadata{jj,3};
                event = panelMetadata{jj,4};
                
                % CHECK IF CELL (EVENT) ADDED TWICE --> REMOVE CELL FROM
                % ANALYSIS
                eventFound = 0;
                for pp = 1:size(table,1)
                    
                    EvalxVar = table2array(table(pp,4));
                    EvalyVar = table2array(table(pp,5));
                    
                    % evaluate if the coordinates belong to the current ROI
                    if EvalxVar >= xLim && EvalxVar <= xLim+xDim   && ...
                            EvalyVar >= yLim && EvalyVar <= yLim+yDim
                        eventFound = eventFound + 1;
                    end
                end
                
                if eventFound >= 2
                    disqualifyCellStatus = 'yes';
                else
                    disqualifyCellStatus = 'no';
                end
                
                eventOrder = {'first' 'second' 'third' 'fourth' 'fifth' 'sixth' 'seventh' 'eighth' 'nineth' 'tenth'};
                for ee = 1:9
                    orderStr = sprintf('event_0%d',ee);
                    if contains(event, orderStr) == 1
                        eventField = eventOrder{ee};
                    end
                end
                
                pw = [ip.dataDir filesep exp filesep ts filesep 'matlabOutput'];
                cd (pw)
                
                position = ts;
                if (position(length(position)) == '#') == 1
                    position = position(1:(length(position)-1));
                end
                filename = char(sprintf('%s%s', position,'_cellDataStructure'));
                load(filename)
                
                switch disqualifyCellStatus
                    case 'yes'
                        cellDataStructure.(exp).(position).(cell).disqualifyCell = 'yes';
                    case 'no'
                        cellDataStructure.(exp).(position).(cell).allEvents.(eventField).validationStatus = 'false event';
                end
                
                
                save(char(filename),'cellDataStructure')
            end
        end
    end
    
    fprintf('\n\nValid events\n');
    
    % Change validation status of non-false events as valid events
    for jj = 1:size(panelMetadata,1)
        exp = panelMetadata{jj,1};
        ts = panelMetadata{jj,2};
        cell = panelMetadata{jj,3};
        event = panelMetadata{jj,4};
        
        eventOrder = {'first' 'second' 'third' 'fourth' 'fifth' 'sixth' 'seventh' 'eighth' 'nineth' 'tenth'};
        for ee = 1:9
            orderStr = sprintf('event_0%d',ee);
            if contains(event, orderStr) == 1
                eventField = eventOrder{ee};
            end
        end
        
        pw = [ip.dataDir filesep exp filesep ts filesep 'matlabOutput'];
        cd (pw)
        
        position = ts;
        if (position(length(position)) == '#') == 1
            position = position(1:(length(position)-1));
        end
        filename = char(sprintf('%s%s', position,'_cellDataStructure'));
        load(filename)
        
        switch cellDataStructure.(exp).(position).(cell).disqualifyCell
            case 'no'
                switch cellDataStructure.(exp).(position).(cell).allEvents.(eventField).validationStatus
                    case 'not validated'
                        cellDataStructure.(exp).(position).(cell).allEvents.(eventField).validationStatus = 'valid event';
                        fprintf('.');
                end
        end
        save(char(filename),'cellDataStructure')
    end
    
elseif ip.noFalseEvents == 1
    
    for jj = 1:size(panelMetadata,1)
        exp = panelMetadata{jj,1};
        ts = panelMetadata{jj,2};
        cell = panelMetadata{jj,3};
        event = panelMetadata{jj,4};
        
        eventOrder = {'first' 'second' 'third' 'fourth' 'fifth' 'sixth' 'seventh' 'eighth' 'nineth' 'tenth'};
        for ee = 1:9
            orderStr = sprintf('event_0%d',ee);
            if contains(event, orderStr) == 1
                eventField = eventOrder{ee};
            end
        end
        
        pw = [ip.dataDir filesep exp filesep ts filesep 'matlabOutput'];
        cd (pw)
        
        position = ts;
        if (position(length(position)) == '#') == 1
            position = position(1:(length(position)-1));
        end
        filename = char(sprintf('%s%s', position,'_cellDataStructure'));
        load(filename)
        
        switch cellDataStructure.(exp).(position).(cell).disqualifyCell
            case 'no'
                switch cellDataStructure.(exp).(position).(cell).allEvents.(eventField).validationStatus
                    case 'not validated'
                        cellDataStructure.(exp).(position).(cell).allEvents.(eventField).validationStatus = 'valid event';
                end
        end
        save(char(filename),'cellDataStructure')
    end
end

fprintf('\n\nValidation status has been updated.\n')
    
end

