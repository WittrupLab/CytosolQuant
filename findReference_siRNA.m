function [ finalFWHMplane ] = findReference_siRNA(cdir)
%% Find maximum mean frame intensity image
cd (cdir)
clear dir
dirData = dir; % Get the data for the current directory
dirData(strncmp({dirData.name}, '.', 1)) = []; % Remove false files
dirData(strncmp({dirData.name}, '~', 1)) = [];
allFiles = {dirData.name}';

nPlanes = numel(allFiles)/2;
imStack = zeros(1024,1024,nPlanes);

addFrame = 1;
for ii = 1:length(allFiles)
    if contains(allFiles(ii),'c2') == 1
        imStack(:,:,addFrame) = imread(char(allFiles(ii)));
        addFrame = addFrame + 1;
    end
end

D = [];
cellSize = size(imStack); %Nr of planes

for g = 1:cellSize(1,3)
    C = sum(imStack(:,:,g)); %Sum pixel values in 1st dimension
    D(g,1) = sum(C); %Sum pixel values in 2nd dimension
end

planeMean = D./(1024*1024); %Calculate plane mean
planeMax = max(planeMean); %Find max plane
FWHM = planeMean(1:end,1)>(0.5*planeMax); %Find planes that meet full width half maximal requirement
FWHM = FWHM.'; %column to row 
FWHMplanes = imStack(:,:,FWHM); %Collect FWHM planes
finalFWHMplane = mean(FWHMplanes,3); %Mean of planes

end

