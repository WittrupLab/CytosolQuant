function [] = runApp(app,task)
%% Use cytoQuant app inputs to initiate processing tasks
ip.homeDir = app.Home.Value;
ip.dataDir = [app.Home.Value filesep 'Data repository'];
ip.fdp = [app.Home.Value filesep 'Processing tools' filesep 'Matlab'];
ip.mdp = [app.Home.Value filesep 'Processing tools' filesep 'Tracking masks' filesep 'masks_170710'];
ip.outputDir = [app.Home.Value filesep 'Data output'];

ip.allExperiments = cell(0);

if app.activateTab.Value == 1
    index = numel(ip.allExperiments)+1;
    ip.allExperiments(index) = {app.experimentID_Tab.Value};
    
    switch app.ChC_Tab.Value
        case 'Not selected'
            ip.allChannels.(ip.allExperiments{index}) = ...
                {app.ChA_Tab.Value app.ChB_Tab.Value};
                  
        otherwise
            ip.allChannels.(ip.allExperiments{index}) = ...
                {app.ChA_Tab.Value app.ChB_Tab.Value app.ChC_Tab.Value};

    end

    ip.nAcquisitions.(ip.allExperiments{index}) = app.TS_Tab.Value;
    
    if app.includeTS_all_Tab.Value == 1
        ip.processAll_TS.(ip.allExperiments{index}) = 1;
        
    elseif app.includeTS_selected_Tab.Value == 1
        ip.processAll_TS.(ip.allExperiments{index}) = 0;
    end
        
    ip.xDim.(ip.allExperiments{index}) = app.xDim_Tab.Value;
    ip.yDim.(ip.allExperiments{index}) = app.yDim_Tab.Value;
    ip.zDim.(ip.allExperiments{index}) = app.zDim_Tab.Value;
    
    ip.defineControls.(ip.allExperiments{index}) = [...
        app.exp1_pos1.Value app.exp1_pos2.Value app.exp1_pos3.Value app.exp1_pos4.Value app.exp1_pos5.Value ...
        app.exp1_pos6.Value app.exp1_pos7.Value app.exp1_pos8.Value app.exp1_pos9.Value app.exp1_pos10.Value]';
end

if app.activateTab_2.Value == 1
    index = numel(ip.allExperiments)+1;
    ip.allExperiments(index) = {app.experimentID_Tab_2.Value};
    
    switch app.ChC_Tab_2.Value
        case 'Not selected'
            ip.allChannels.(ip.allExperiments{index}) = ...
                {app.ChA_Tab_2.Value app.ChB_Tab_2.Value};
            
        otherwise
            ip.allChannels.(ip.allExperiments{index}) = ...
                {app.ChA_Tab_2.Value app.ChB_Tab_2.Value app.ChC_Tab_2.Value};

    end

    ip.nAcquisitions.(ip.allExperiments{index}) = app.TS_Tab_2.Value;
    
    if app.includeTS_all_Tab_2.Value == 1
        ip.processAll_TS.(ip.allExperiments{index}) = 1;
        
    elseif app.includeTS_selected_Tab_2.Value == 1
        ip.processAll_TS.(ip.allExperiments{index}) = 0;
    end
        
    ip.xDim.(ip.allExperiments{index}) = app.xDim_Tab_2.Value;
    ip.yDim.(ip.allExperiments{index}) = app.yDim_Tab_2.Value;
    ip.zDim.(ip.allExperiments{index}) = app.zDim_Tab_2.Value;
    
    ip.defineControls.(ip.allExperiments{index}) = [...
        app.exp2_pos1.Value app.exp2_pos2.Value app.exp2_pos3.Value app.exp2_pos4.Value app.exp2_pos5.Value ...
        app.exp2_pos6.Value app.exp2_pos7.Value app.exp2_pos8.Value app.exp2_pos9.Value app.exp2_pos10.Value]';
end

if app.activateTab_3.Value == 1
    index = numel(ip.allExperiments)+1;
    ip.allExperiments(index) = {app.experimentID_Tab_3.Value};
    
    switch app.ChC_Tab_3.Value
        case 'Not selected'
            ip.allChannels.(ip.allExperiments{index}) = ...
                {app.ChA_Tab_3.Value app.ChB_Tab_3.Value}; 
            
        otherwise
            ip.allChannels.(ip.allExperiments{index}) = ...
                {app.ChA_Tab_3.Value app.ChB_Tab_3.Value app.ChC_Tab_3.Value};
    end

    ip.nAcquisitions.(ip.allExperiments{index}) = app.TS_Tab_3.Value;
    
    if app.includeTS_all_Tab_3.Value == 1
        ip.processAll_TS.(ip.allExperiments{index}) = 1;
        
    elseif app.includeTS_selected_Tab_3.Value == 1
        ip.processAll_TS.(ip.allExperiments{index}) = 0;
    end
        
    ip.xDim.(ip.allExperiments{index}) = app.xDim_Tab_3.Value;
    ip.yDim.(ip.allExperiments{index}) = app.yDim_Tab_3.Value;
    ip.zDim.(ip.allExperiments{index}) = app.zDim_Tab_3.Value;
    
    ip.defineControls.(ip.allExperiments{index}) = [...
        app.exp3_pos1.Value app.exp3_pos2.Value app.exp3_pos3.Value app.exp3_pos4.Value app.exp3_pos5.Value ...
        app.exp3_pos6.Value app.exp3_pos7.Value app.exp3_pos8.Value app.exp3_pos9.Value app.exp3_pos10.Value]';
    
end

if app.activateTab_4.Value == 1
    index = numel(ip.allExperiments)+1;
    ip.allExperiments(index) = {app.experimentID_Tab_4.Value};
    
    switch app.ChC_Tab_4.Value
        case 'Not selected'
            ip.allChannels.(ip.allExperiments{index}) = ...
                {app.ChA_Tab_4.Value app.ChB_Tab_4.Value};   
            
        otherwise
            ip.allChannels.(ip.allExperiments{index}) = ...
                {app.ChA_Tab_4.Value app.ChB_Tab_4.Value app.ChC_Tab_4.Value};
    end

    ip.nAcquisitions.(ip.allExperiments{index}) = app.TS_Tab_4.Value;
    
    if app.includeTS_all_Tab_4.Value == 1
        ip.processAll_TS.(ip.allExperiments{index}) = 1;
        
    elseif app.includeTS_selected_Tab_4.Value == 1
        ip.processAll_TS.(ip.allExperiments{index}) = 0;
    end
        
    ip.xDim.(ip.allExperiments{index}) = app.xDim_Tab_4.Value;
    ip.yDim.(ip.allExperiments{index}) = app.yDim_Tab_4.Value;
    ip.zDim.(ip.allExperiments{index}) = app.zDim_Tab_4.Value;
    
    ip.defineControls.(ip.allExperiments{index}) = [...
        app.exp4_pos1.Value app.exp4_pos2.Value app.exp4_pos3.Value app.exp4_pos4.Value app.exp4_pos5.Value ...
        app.exp4_pos6.Value app.exp4_pos7.Value app.exp4_pos8.Value app.exp4_pos9.Value app.exp4_pos10.Value]';
    
end

if app.activateTab_5.Value == 1
    index = numel(ip.allExperiments)+1;
    ip.allExperiments(index) = {app.experimentID_Tab_5.Value};
    
    switch app.ChC_Tab_5.Value
        case 'Not selected'
            ip.allChannels.(ip.allExperiments{index}) = ...
                {app.ChA_Tab_5.Value app.ChB_Tab_5.Value}; 
            
        otherwise
            ip.allChannels.(ip.allExperiments{index}) = ...
                {app.ChA_Tab_5.Value app.ChB_Tab_5.Value app.ChC_Tab_5.Value};
    end

    ip.nAcquisitions.(ip.allExperiments{index}) = app.TS_Tab_5.Value;
    
    if app.includeTS_all_Tab_5.Value == 1
        ip.processAll_TS.(ip.allExperiments{index}) = 1;
        
    elseif app.includeTS_selected_Tab_5.Value == 1
        ip.processAll_TS.(ip.allExperiments{index}) = 0;
    end
        
    ip.xDim.(ip.allExperiments{index}) = app.xDim_Tab_5.Value;
    ip.yDim.(ip.allExperiments{index}) = app.yDim_Tab_5.Value;
    ip.zDim.(ip.allExperiments{index}) = app.zDim_Tab_5.Value;
    
    ip.defineControls.(ip.allExperiments{index}) = [...
        app.exp5_pos1.Value app.exp5_pos2.Value app.exp5_pos3.Value app.exp5_pos4.Value app.exp5_pos5.Value ...
        app.exp5_pos6.Value app.exp5_pos7.Value app.exp5_pos8.Value app.exp5_pos9.Value app.exp5_pos10.Value]';
    
end

if app.activateTab_6.Value == 1
    index = numel(ip.allExperiments)+1;
    ip.allExperiments(index) = {app.experimentID_Tab_6.Value};
    
    switch app.ChC_Tab_6.Value
        case 'Not selected'
            ip.allChannels.(ip.allExperiments{index}) = ...
                {app.ChA_Tab_6.Value app.ChB_Tab_6.Value};
            
        otherwise
            ip.allChannels.(ip.allExperiments{index}) = ...
                {app.ChA_Tab_6.Value app.ChB_Tab_6.Value app.ChC_Tab_6.Value};
    end

    ip.nAcquisitions.(ip.allExperiments{index}) = app.TS_Tab_6.Value;
    
    if app.includeTS_all_Tab_6.Value == 1
        ip.processAll_TS.(ip.allExperiments{index}) = 1;
        
    elseif app.includeTS_selected_Tab_6.Value == 1
        ip.processAll_TS.(ip.allExperiments{index}) = 0;
    end
        
    ip.xDim.(ip.allExperiments{index}) = app.xDim_Tab_6.Value;
    ip.yDim.(ip.allExperiments{index}) = app.yDim_Tab_6.Value;
    ip.zDim.(ip.allExperiments{index}) = app.zDim_Tab_6.Value;
    
    ip.defineControls.(ip.allExperiments{index}) = [...
        app.exp6_pos1.Value app.exp6_pos2.Value app.exp6_pos3.Value app.exp6_pos4.Value app.exp6_pos5.Value ...
        app.exp6_pos6.Value app.exp6_pos7.Value app.exp6_pos8.Value app.exp6_pos9.Value app.exp6_pos10.Value]';
    
end

if app.activateTab_7.Value == 1
    index = numel(ip.allExperiments)+1;
    ip.allExperiments(index) = {app.experimentID_Tab_7.Value};
    
    switch app.ChC_Tab_7.Value
        case 'Not selected'
            ip.allChannels.(ip.allExperiments{index}) = ...
                {app.ChA_Tab_7.Value app.ChB_Tab_7.Value};
            
        otherwise
            ip.allChannels.(ip.allExperiments{index}) = ...
                {app.ChA_Tab_7.Value app.ChB_Tab_7.Value app.ChC_Tab_7.Value};
    end

    ip.nAcquisitions.(ip.allExperiments{index}) = app.TS_Tab_7.Value;
    
    if app.includeTS_all_Tab_7.Value == 1
        ip.processAll_TS.(ip.allExperiments{index}) = 1;
        
    elseif app.includeTS_selected_Tab_7.Value == 1
        ip.processAll_TS.(ip.allExperiments{index}) = 0;
    end
        
    ip.xDim.(ip.allExperiments{index}) = app.xDim_Tab_7.Value;
    ip.yDim.(ip.allExperiments{index}) = app.yDim_Tab_7.Value;
    ip.zDim.(ip.allExperiments{index}) = app.zDim_Tab_7.Value;
    
    ip.defineControls.(ip.allExperiments{index}) = [...
        app.exp7_pos1.Value app.exp7_pos2.Value app.exp7_pos3.Value app.exp7_pos4.Value app.exp7_pos5.Value ...
        app.exp7_pos6.Value app.exp7_pos7.Value app.exp7_pos8.Value app.exp7_pos9.Value app.exp7_pos10.Value]';
end

if app.activateTab_8.Value == 1
    index = numel(ip.allExperiments)+1;
    ip.allExperiments(index) = {app.experimentID_Tab_8.Value};
    
    switch app.ChC_Tab_8.Value
        case 'Not selected'
            ip.allChannels.(ip.allExperiments{index}) = ...
                {app.ChA_Tab_8.Value app.ChB_Tab_8.Value}; 
            
        otherwise
            ip.allChannels.(ip.allExperiments{index}) = ...
                {app.ChA_Tab_8.Value app.ChB_Tab_8.Value app.ChC_Tab_8.Value};
    end


    
    ip.nAcquisitions.(ip.allExperiments{index}) = app.TS_Tab_8.Value;
    
    if app.includeTS_all_Tab_8.Value == 1
        ip.processAll_TS.(ip.allExperiments{index}) = 1;
        
    elseif app.includeTS_selected_Tab_8.Value == 1
        ip.processAll_TS.(ip.allExperiments{index}) = 0;
    end
        
    ip.xDim.(ip.allExperiments{index}) = app.xDim_Tab_8.Value;
    ip.yDim.(ip.allExperiments{index}) = app.yDim_Tab_8.Value;
    ip.zDim.(ip.allExperiments{index}) = app.zDim_Tab_8.Value;
    
    ip.defineControls.(ip.allExperiments{index}) = [...
        app.exp8_pos1.Value app.exp8_pos2.Value app.exp8_pos3.Value app.exp8_pos4.Value app.exp8_pos5.Value ...
        app.exp8_pos6.Value app.exp8_pos7.Value app.exp8_pos8.Value app.exp8_pos9.Value app.exp8_pos10.Value]';
end

if app.activateTab_9.Value == 1
    index = numel(ip.allExperiments)+1;
    ip.allExperiments(index) = {app.experimentID_Tab_9.Value};
    
    switch app.ChC_Tab_9.Value
        case 'Not selected'
            ip.allChannels.(ip.allExperiments{index}) = ...
                {app.ChA_Tab_9.Value app.ChB_Tab_9.Value};    
            
        otherwise
            ip.allChannels.(ip.allExperiments{index}) = ...
                {app.ChA_Tab_9.Value app.ChB_Tab_9.Value app.ChC_Tab_9.Value};
    end
 
    ip.nAcquisitions.(ip.allExperiments{index}) = app.TS_Tab_9.Value;
    
    if app.includeTS_all_Tab_9.Value == 1
        ip.processAll_TS.(ip.allExperiments{index}) = 1;
        
    elseif app.includeTS_selected_Tab_9.Value == 1
        ip.processAll_TS.(ip.allExperiments{index}) = 0;
    end
        
    ip.xDim.(ip.allExperiments{index}) = app.xDim_Tab_9.Value;
    ip.yDim.(ip.allExperiments{index}) = app.yDim_Tab_9.Value;
    ip.zDim.(ip.allExperiments{index}) = app.zDim_Tab_9.Value;
    
    ip.defineControls.(ip.allExperiments{index}) = [...
        app.exp9_pos1.Value app.exp9_pos2.Value app.exp9_pos3.Value app.exp9_pos4.Value app.exp9_pos5.Value ...
        app.exp9_pos6.Value app.exp9_pos7.Value app.exp9_pos8.Value app.exp9_pos9.Value app.exp9_pos10.Value]';
end

if app.activateTab_10.Value == 1
    index = numel(ip.allExperiments)+1;
    ip.allExperiments(index) = {app.experimentID_Tab_10.Value};
    
    switch app.ChC_Tab_10.Value
        case 'Not selected'
            ip.allChannels.(ip.allExperiments{index}) = ...
                {app.ChA_Tab_10.Value app.ChB_Tab_10.Value};   
            
        otherwise
            ip.allChannels.(ip.allExperiments{index}) = ...
                {app.ChA_Tab_10.Value app.ChB_Tab_10.Value app.ChC_Tab_10.Value};
    end

    ip.nAcquisitions.(ip.allExperiments{index}) = app.TS_Tab_10.Value;
    
    if app.includeTS_all_Tab_10.Value == 1
        ip.processAll_TS.(ip.allExperiments{index}) = 1;
        
    elseif app.includeTS_selected_Tab_10.Value == 1
        ip.processAll_TS.(ip.allExperiments{index}) = 0;
    end
        
    ip.xDim.(ip.allExperiments{index}) = app.xDim_Tab_10.Value;
    ip.yDim.(ip.allExperiments{index}) = app.yDim_Tab_10.Value;
    ip.zDim.(ip.allExperiments{index}) = app.zDim_Tab_10.Value;
    
    ip.defineControls.(ip.allExperiments{index}) = [...
        app.exp10_pos1.Value app.exp10_pos2.Value app.exp10_pos3.Value app.exp10_pos4.Value app.exp10_pos5.Value ...
        app.exp10_pos6.Value app.exp10_pos7.Value app.exp10_pos8.Value app.exp10_pos9.Value app.exp10_pos10.Value]';
end

switch task
    case 'Calculate bleach correction'
        calcBleachCorrection(ip)
    
    case 'Detect release events'
        
        [mitosisEvents] = analyzeMultiPosKnockdownPrimtoSec(ip);
        assignin('base', 'mitosisEvents', mitosisEvents);
    
    case 'Create folders'
        cd (ip.fdp)
        createFolders(ip)

    case 'Primary processing'
        pause on
        delayStart = 0; %app.DelayStartSlider.Value;
        pause(60*60*delayStart)
             
        cd (ip.fdp)
        ip.skipDecon = 1;
        primaryProcessing(ip);
        
    case 'Secondary processing'
        cd (ip.fdp)
        secondaryProcessing(ip);
        
    case 'Create ROIs'
        ip.extROI = app.ExtendROI.Value;
        ip.extTimePre = app.TimesPre.Value;
        ip.extTimePost = app.TimesPost.Value;
        ip.executionMode = 'New';
        ip.fileAccessMode = 'NAS';

        if app.roiTS_all.Value == 1
           ip.roiAll_TS = 1;
        elseif app.roiTS_selected.Value == 1
           ip.roiAll_TS = 0;
        end
        
        if app.roiEvents_all.Value == 1
           ip.roiAll_events = 1;
        elseif app.roiEvents_selected.Value == 1
           ip.roiAll_events = 0;
        end

       ip.skipCelldata = 1;
       ip.skipImageBackground = 1;
        
        cd (ip.fdp)
        createROIs(ip)
                 
    case 'Define controls'
        defineControls(ip);
        
    case 'Collect tracking data'
        if app.collectTS_all.Value == 1
             ip.processAll_TS =  1;
        elseif app.collectTS_selected.Value == 1
             ip.processAll_TS =  0;
        end
        
        if app.collectEvents_all.Value == 1
              ip.processAll_events = 1;
        elseif app.collectEvents_selected.Value == 1
              ip.processAll_events = 0;
        end

        ip.collect.controlCells = app.collect_controlCells.Value;
        ip.collect.nonEventCells = app.collect_nonEvenCells.Value;
        ip.collect.mitosisCells_nonEvenCells_timeShifted = app.collect_mitosisCells_nonEvenCells_timeShifted.Value;
        ip.collect.apoptoticNonEventCells =  app.collect_apoptoticNonEventCells;
        ip.collect.eventCells_timeShifted_primary =  app.collect_eventCells_timeShifted_primary.Value;
        ip.collect.eventCells_timeShifted_secondary =  app.collect_eventCells_timeShifted_secondary.Value;
        ip.collect.eventCells_timeShifted_tertiary = app.collect_eventCells_timeShifted_tertiary.Value;
        ip.collect.eventCells_timeShifted_quatary = app.collect_eventCells_timeShifted_quatery.Value;
        ip.collect.mitosisCells_evenCells_timeShifted = app.collect_mitosisCells_evenCells_timeShifted.Value;
        ip.collect.nCellsWithEvents = app.collect_nCellsWithEvents.Value;
        ip.collect.timediff_event_mitosis = app.collect_timediff_event_mitosis.Value;
        ip.collect.apoptoticEventCells = app.collect_apoptoticEventCells.Value;
        ip.collect.normGFP = app.collect_normGFP.Value;

        cd (ip.fdp)
        [collectedData] = collectData(ip);
        assignin('base', 'collectedData', collectedData)
        
    case 'Create event panel'
        ip.mode = 'mipUntracked';
     
        ip.extPlanes = app.ExtendPlanes.Value;
        ip.spacersx = app.xSpacers.Value;
        ip.spacersy = app.ySpacers.Value;
        ip.frameSizeY = app.yPanelDim.Value;
        ip.frameSizeX = app.xPanelDim.Value;
        
        if app.eventPanelTS_all.Value == 1
            ip.plotAll_TS = 1;
        elseif app.eventPanelTS_selected.Value == 1
            ip.plotAll_TS = 0;
        end
        
        if app.eventPanelCells_all.Value == 1
            ip.plotAll_cells =  1;
        elseif app.eventPanelCells_selected.Value == 1
            ip.plotAll_cells =  0;
        end
        
        if app.eventPanelEvents_all.Value == 1
            ip.plotAll_events =   1;
        elseif app.eventPanelEvents_selected.Value == 1
            ip.plotAll_events =  0;
        end
        
        ip.writepath = app.exportFolder.Value;
        ip.panelName = app.panelName.Value;
        
        if app.notValidatedEvents.Value == 1
            ip.includeEventStatus = 'not validated';
        elseif app.validatedEvents.Value == 1
            ip.includeEventStatus = 'valid events';
        elseif app.falseEvents.Value == 1
            ip.includeEventStatus = 'false events';
        else
            ip.includeEventStatus = 'all events';
        end
        
        eventOrder = {'first' 'second' 'third' 'fourth' 'fifth' 'sixth' 'seventh' 'eighth' 'nineth' 'tenth'};
        for ee = 1:9
            selectOrder = app.includeEventOrder.Value;
            loopIndex = sprintf('%d',ee);
            if contains(selectOrder, loopIndex)
                ip.includeEventOrder = eventOrder{ee};
            end
        end
        
        cd (ip.fdp)
        eventPanel(ip);

    case 'Export panel'
        ip.writepath = app.exportFolder.Value;
        ip.panelName = app.panelName.Value;
        cd (ip.fdp)
        exportEventPanel(panel,panelMetadata,writepath,panelName)

    case 'Add validation data'
        readpath = app.panelFolderPath.Value;
        eventOrder = {'first' 'second' 'third' 'fourth' 'fifth' 'sixth' 'seventh' 'eighth' 'nineth' 'tenth'};
        for ee = 1:9
            selectOrder = app.includeEventOrder.Value;
            loopIndex = sprintf('%d',ee);
            if contains(selectOrder, loopIndex)
                ip.includeEventOrder = eventOrder{ee};
            end
        end
        
        if app.noFalseEvents.Value == 0
            ip.noFalseEvents = 0;
        else
            ip.noFalseEvents = 1;
        end

        if app.notValidatedEvents.Value == 1
            ip.includeEventStatus = 'not validated';
        elseif app.validatedEvents.Value == 1
            ip.includeEventStatus = 'valid events';
        elseif app.falseEvents.Value == 1
            ip.includeEventStatus = 'false events';
        else
            ip.includeEventStatus = 'all events';
        end
        
        cd (ip.fdp)
        addValidationData(ip,readpath)
end

end

