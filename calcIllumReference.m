function [illumReference] = calcIllumReference(ip,cdir)
%% Calculate illumination reference image from siRNA and control images
cd (cdir)
cd 'referenceMeasure/rawTiffs/P11000nMAF6473343v2'
cdir = pwd;
cd (ip.fdp)

siRNA_ref1 = findReference_siRNA(cdir);

cd (cdir)
cd ../
cd 'P21000nMAF6473343v2'
cdir = pwd;
cd (ip.fdp)

siRNA_ref2 = findReference_siRNA(cdir);


cd (cdir)
cd ../
cd 'P31000nMAF6473343v2'
cdir = pwd;
cd (ip.fdp)

siRNA_ref3 = findReference_siRNA(cdir);


cd (cdir)
cd ../
cd 'P1CMBCTRL'
cdir = pwd;
cd (ip.fdp)

CTRL_ref1 = findReference_CTRL(cdir);


cd (cdir)
cd ../
cd 'P2CMBCTRL'
cdir = pwd;
cd (ip.fdp)

CTRL_ref2 = findReference_CTRL(cdir);


cd (cdir)
cd ../
cd 'P3CMBCTRL'
cdir = pwd;
cd (ip.fdp)

CTRL_ref3 = findReference_CTRL(cdir);


cd (cdir)
cd ../../..


siRNA_stack = siRNA_ref1;
siRNA_stack(:,:,2) = siRNA_ref2;
siRNA_stack(:,:,3) = siRNA_ref3;
siRNA_mean = mean(siRNA_stack,3);

CTRL_stack = CTRL_ref1;
CTRL_stack(:,:,2) = CTRL_ref2;
CTRL_stack(:,:,3) = CTRL_ref3;
CTRL_mean = mean(CTRL_stack,3);

illumReference = (siRNA_mean - CTRL_mean) / 1000;

meanFilter = fspecial('disk',20);
illumReference=filter2(meanFilter, illumReference);

end

