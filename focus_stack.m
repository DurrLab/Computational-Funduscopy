% This script shows an example of refocusing a raw diffuser image
% Deconvolves 1 raw images with a range of PSFs at different depths.

% Import data
load("data/calibrationparams.mat","P","RefractiveErrors","deltaZ");
load("data/PSFinfo.mat","pinhole_um","x_samples","y_samples","z_samples");

% Convert image indices to zstack index.
zstack = compute_indices(x_samples, y_samples, z_samples);
zlen = length(z_samples);
psfdir = "data/PSFs/";
outputdir = "data/Output";

% Raw diffuser image
im = imread("data/RawImages/ModelEye0D.tif");

f = waitbar(0,'Loading data','Name','Reconstructing focal stack');
for k = 300:301 % Loop through indices you want to process
    waitbar(double(k)/zlen, f, sprintf('Frame: %d/%d',k,zlen))
    idx = zstack(k);
    imdir = fullfile(psfdir,strcat(string(idx),'.tif'));
    psf = imread(imdir);
    reconim = recon(im, psf, [0.004, 0.002, 0.002]);
    imwrite(reconim,fullfile(outputdir, strcat(string(k),'.tif')),'tif');
end
delete(f);

%% Functions
% Helper functions
function zstack = compute_indices(x_samples, y_samples, z_samples)
xlen = length(x_samples);
ylen = length(y_samples);
zlen = length(z_samples);
idxarray = zeros(xlen*ylen*zlen, 3);
k = 1;
for x = 1:length(x_samples)
    for y = 1:length(y_samples)
        for z = 1:length(z_samples)
            idxarray(k, :) = [x, y, z];
            k = k + 1;
        end
    end
end
zstack = find(idxarray(:,1) == 1 & idxarray(:,2) == 1);
end
