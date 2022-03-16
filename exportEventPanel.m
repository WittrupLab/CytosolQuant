function [ ] = exportEventPanel(ip, panel, panelMetadata)
%% Save event panel
cd (char(ip.writepath))
destination = 'eventPanels';
if ~exist(destination,'dir')
    mkdir(destination);
end
cd (destination)

dateparam = clock;
destination = char(sprintf('%s %d.%d.%d %d.%d', ip.panelName, dateparam(1),dateparam(2),dateparam(3),dateparam(4),dateparam(5)));
if ~exist(destination,'dir')
    mkdir(destination);
end
cd (destination)

mkdir 'Panel files'
cd 'Panel files'

syntax = 'eventPanel_c%s_t%s.tif';

for iChannel = 1:size(panel,5)
    for iTime = 1:size(panel,4)

        stringC = '%d';
        if iTime < 10
            stringT = '0%d';
        elseif iTime < 100
            stringT = '%d';
        end
        filesyntax = sprintf(char(syntax),stringC, stringT);
        filename = sprintf(filesyntax, iChannel, iTime);
            
        for iPlane = 1:size(panel,3)
            
            data = panel(:,:,iPlane,iTime,iChannel);
            data = uint16(data);
            imwrite(data, char(filename), 'WriteMode', 'append');
        end
    end
end

cd ../

save('panelMetadata.mat', 'panelMetadata');
end

