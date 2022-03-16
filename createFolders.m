function [] = createFolders(ip)
%% Create initial experiment folder architecture   
nExperiments = numel(ip.allExperiments);
cd (char(ip.dataDir));

fprintf('\nCreating folders...\n')

% Create folder architecture in data directory
for ee = 1:nExperiments
    mkdir (char(ip.allExperiments(ee)))
    cd (char(ip.allExperiments(ee)))
    nDirs = ip.nAcquisitions.(char(ip.allExperiments(ee)));
    allChannels = ip.allChannels.(char(ip.allExperiments(ee)));
    
    for tsDirs = 1:nDirs
        if tsDirs < 10
            string1 = '0%d';
        else
            string1 = '%d';
        end
        syntax = sprintf('%s_TS%s', char(ip.allExperiments(ee)), string1);
        
        foldername = sprintf(char(syntax),tsDirs);
        mkdir (char(foldername));
        cd (char(foldername))
        
        mkdir ('zenRaw');
        mkdir ('rawTiffs');
        mkdir ('denoisedNuc');
        mkdir ('denoisedGFP');
        mkdir ('cellprofilerOutput');
        
        
        cd ('cellprofilerOutput');
        mkdir ('segNuc');
        mkdir ('segCell');
        mkdir ('segRelabled');
        mkdir ('dataSheet');
        
        cd ('dataSheet');
        mkdir ('identifyNuc');
        mkdir ('cellFullsiRNA');
        mkdir ('cellMaskedsiRNA');
        mkdir ('backgroundIntensity');
        
        cd ../..
        
        mkdir ('matlabOutput');
        
        if tsDirs == 1
            mkdir ('cellProfilerPipelines')
            mkdir ('referenceMeasure');
            cd ('referenceMeasure');
            mkdir ('zenRaw');
            mkdir ('rawTiffs');
            mkdir ('cellprofilerOutput');
            cd ('cellprofilerOutput');
            mkdir ('1000nM');
            mkdir ('0nM');
            cd ../..
        end
        
        mkdir ('readMetadata');
        mkdir ('eventMetadata');
        mkdir ('eventROIs');
        mkdir ('readImages');
        mkdir ('falseEvents')
        
        cd ../
    end
    cd ../
end

fprintf('\nFolders have been created.\n')

end

