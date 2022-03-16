function [] = defineControls(ip)
%% Define control positions
fprintf('\n\nDefining control positions...\n');
cd (ip.fdp)
ip.listAcquisitions = listAcquisitions(ip);

for iAcquisitions = 1:size(ip.listAcquisitions,1)
    ip.iAcquisitions = iAcquisitions;
    experimentFolder = [ip.dataDir filesep ip.listAcquisitions{ip.iAcquisitions,1}];
    fprintf('Processing %s\n',ip.listAcquisitions{iAcquisitions,2});
    position = ip.listAcquisitions{ip.iAcquisitions,2};
    if (position(length(position)) == '#') == 1
        position = position(1:(length(position)-1));
    end
    
    P = {'TS01' 'TS02' 'TS03' 'TS04' 'TS05' ...
         'TS06' 'TS07' 'TS08' 'TS09' 'TS10' };
    controlPositionLogic = ip.defineControls.(ip.listAcquisitions{ip.iAcquisitions,1})';
    controlPositionLogicTable = array2table(controlPositionLogic);
    controlPositionLogicTable.Properties.VariableNames = P;
    
    folder = [experimentFolder filesep ip.listAcquisitions{ip.iAcquisitions,2} filesep 'matlabOutput'];
    cd (folder)
    save(sprintf('%s%s', position, '_defineControls'),'controlPositionLogicTable');
end
cd (ip.fdp)

fprintf('\nControl positions status have been added.\n');
end

