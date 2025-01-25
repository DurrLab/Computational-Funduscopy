% Script to analyze the resolution of the system using PSF MTFs.
%   The maximum theoretical resolution is defined here to be the highest 
%   spatial frequency in the radially averaged MTF that is above the noise
%   floor of the camera.

% Dependencies:      - RayMat class (included in repository)
%                    - radialavg.zip 
%                       (David Fischer (2025). radialavg.zip
%                       (https://www.mathworks.com/matlabcentral/fileexchange/46468-radialavg-zip)
%                       MATLAB Central File Exchange. 
%                     - shadedErrorBar
%                          Rob Campbell (2025). raacampbell/shadedErrorBar 
%                          (https://github.com/raacampbell/shadedErrorBar), 
%                          GitHub. 
%
% Variables:
%   P                -   Optical power of the eye used (diopters). We used 39 diopter
%                        (1 inch focal length) lens. A human eye is 60
%                        diopters. This causes a simple magnification
%                        difference. Using 60 diopters in the code takes
%                        care of this.
%   RefractiveErrors -   Defocus error (in diopters) corresponding to each
%                        image in PSFs
%   deltaZ           -   Defocus distance of PSF (in mm) corresponding to
%                        each image in PSFs
%   pinhole_um       -   Pinhole size (microns)
%   x_samples        -   Horizontal displacements of point source in PSFs
%   y_samples        -   Vertical displacements of point source in PSFs
%   z_samples        -   Axial displacements of point source in PSFs

% Author: 
%    Corey Simmerer

% Download PSFs into .data/PSFs
% Calibrationparams.mat and PSFinfo.mat are in repository.

% Import data
load("data/calibrationparams.mat","P","RefractiveErrors","deltaZ");
load("data/PSFinfo.mat","pinhole_um","x_samples","y_samples","z_samples");

% Convert image indices to zstack index.
zstack = compute_indices(x_samples, y_samples, z_samples);
psfdir = "data/PSFs/";

%% Variables to store results in
tested_k = 1:801; % Indices of the PSFs to test
tested_defocus_errors = RefractiveErrors(tested_k); % Defocus errors for each k
maximum_theoretical_resolutions = zeros(1,length(tested_k)); % Maximum theoretical resolution
error_bars = zeros(1,length(tested_k));
pixel_sizes = zeros(1,length(tested_k));

 
w = waitbar(0,"");
for k_idx = 1:length(tested_k) % specific PSF indices to test
    k = tested_k(k_idx);
    % Defocus error of this PSF index, k
    waitbar(k_idx./length(tested_k),w,string(k_idx));
    eps = defocusError(deltaZ(k),P); % P is the model eye diopters
    dz  = dz_from_eps(eps, 60); % Get the eye-equivalent dz from the defocus (Assume a 60D eye)
    dz = dz*1000; % Convert dz to mm
    f_eye = (1/60)*1000; % focal length of the 60D human eye

    %% Simulate system to get the pixel size on the object
    [~, ~, R, ~] = compute_mag_res(dz, f_eye); % Get pixel size of system (mm)
    pixel_sizes(k_idx) = R;
    %% Get the PSF for this index
    idx = zstack(k);
    imdir = fullfile(psfdir,strcat(string(idx),'.tif'));
    psf = imread(imdir);
    %% Compute the MTF
    otf = fftshift(fft2(psf));
    mtf = abs(otf);
    norm_mtf = mtf./max(mtf(:));
    %% Radially average the MTF
    [avg, tics] = radialavg(norm_mtf, 1024);
    %% Compute the frequency steps in the FFT
     T = R./1000; %30e-6; % 30 micron pixel pitch (unit m)
    fs = 1/(T*1000); % 1 samples per mm
    L = 2048; % Length of signal
    freq_step = fs/L;
    freq_sampling = linspace(-fs/2,fs/2, L);
    radial_freqs = freq_step.*(1:length(tics)); % List of the frequencies
    %% Find the max theoretical resolution
    max_contrast = 1/26667; % Funciton of camera (FWC/noise)
    idx_maxres = find(avg<max_contrast);
    idx_maxres = min(idx_maxres(idx_maxres>1));
    maxres = radial_freqs(idx_maxres);
    
    % Now, we find a linear fit about this point where it crosses the max
    % contrast.
    local_idx = idx_maxres-18:1:idx_maxres+18;
    local_freqs = radial_freqs(local_idx);
    local_mtf = avg(local_idx);
    Lfit = polyfit(local_freqs, local_mtf,1); %P(1) is m, P(2) is b
    m = Lfit(1);

    % Find the standard deviation of the residuals for error bar
    residuals = local_mtf - (Lfit(1)*local_freqs + Lfit(2));
    stddev = std(residuals);

    % Theoretical resolution is where line of best fit crosses 
    maximum_theoretical_resolution = (max_contrast-Lfit(2))/m;
    maximum_theoretical_resolutions(k_idx)=maximum_theoretical_resolution;

    % Error bar length =  std / sin(atan(m))
    error_bar = stddev / sin(atan(m));
    error_bars(k_idx) = error_bar;

end

%% Create plots
% Smoothed version
figure(1); clf;
yyaxis right; hold on;
shadedErrorBar(tested_defocus_errors, smoothdata(maximum_theoretical_resolutions,"gaussian",25),smoothdata(2*error_bars,"movmean",10),lineProps={"Color",[0.8500, 0.3250, 0.0980],"LineWidth",2});
xlabel('Defocus error (D)')
ylabel('Max. resolution (lp/mm)')
xlim([-12 12])
yyaxis left;
plot(tested_defocus_errors, pixel_sizes*1000, LineWidth=2);
ylabel('Size of pixel on fundus (μm)')

%% -----Functions-----
%% Create system model
function system = model_system(dz, f_eye, d_DS)
%   model_system Creates a RayMat system to model the funduscope
%   inputs: 
%       dz     - displacement of point source (mm)
%       f_eye  - focal length of the eye in the system (mm)
%       d_DS   - diffuser-sensor distance (mm)
%   outputs:
%       system - RayMat class representing the system

% Constant parameters (all in mm):
f_L1 = 46.2; 
f_L2 = 75; 
f_L3 = 100; 
dPL1 = 70.85; 
dL1HM = 132.8;
dHML2 = 68.75;
dL2L3 = 96.75;

% 1: Propagation from retina to cornea
R2C = RayMat.freespace(f_eye-dz);
% 2: Cornea gives 60D of refraction
C = RayMat.thinlens(f_eye);
% 3: Prop from pupil to L1
C2L1 = RayMat.freespace(dPL1);
% 4: L1
L1 = RayMat.thinlens(f_L1);
% 5: L1 to HM
L12HM = RayMat.freespace(dL1HM);
% 6: HM to L2
HM2L2 = RayMat.freespace(dHML2);
% 7: L2
L2 = RayMat.thinlens(f_L2);
% 8: L2 to L3
L22L3 = RayMat.freespace(dL2L3);
% 9: L3
L3 = RayMat.thinlens(f_L3);
% 10: L3 to Rpp
L32Rpp = RayMat.freespace(32.52);
% 11: Rpp to Cpp
Rpp2Cpp = RayMat.freespace(79.2);
% 12: Cpp to S
Cpp2S = RayMat.freespace(d_DS);

system_elements = {R2C, C, C2L1, L1, L12HM, HM2L2, L2, L22L3, L3, L32Rpp, Rpp2Cpp, Cpp2S};
system = RayMat(system_elements);
end

function [d_L3_Fpp, M_Lat, R, res] = compute_mag_res(dz, f_eye)
    % compute_mag_res computes the magnification and smallest resolvable
    % feature (pixel size in object space) of the system with defocus
    % caused by dz.
    % Inputs:  dz       - defocus distance (mm)
    %          f_eye    - focal length of the eye used (mm)
    % Outputs: d_L3_Fpp - Distance from lens 3 to F'', the final fundus
    %                     conjugate plane.
    %          M_Lat    - Lateral resolution of system
    %          R        - Smallest resolvable feature using pixel size in 
    %                     object space (mm)
    %          res      - Nyquist limit given pixel size R (lp/mm)

    d_DS = 11; % mm, constant
    fullSystem = model_system(dz, f_eye, d_DS);
    C2L3System = RayMat({fullSystem.elements{2:9}}); % Subsystem from Cornea to L3
    % This is the [A B C D] matrix from Eq (5) in paper.
    a = C2L3System.rayMatrix(1,1);
    b = C2L3System.rayMatrix(1,2);
    c = C2L3System.rayMatrix(2,1);
    d = C2L3System.rayMatrix(2,2);
    d_L3_Fpp = -(a.*f_eye-a.*dz+b)./(c.*f_eye-c.*dz+d); % distance from L3 to F'' Eq (6) in paper
    M_Lat = a + c.*d_L3_Fpp; % Lateral magnifcation. Eq (7) in paper.
    d_L3_D = 111.72; % fixed diffuser-L3 distance.
    d_Fpp_D = d_L3_D - d_L3_Fpp;
    pixel = 6.4/1000; % 6.4 um. (unit: mm) Pixel pitch of PCO 4.2bi uv
    R = pixel * (1/M_Lat) * (d_Fpp_D/d_DS); % Smallest resolvable feautre (pixel size in object space) in mm.
    res = 1/(2*R);
end

%% Useful formulas
function epsilon = defocusError(dz, P)
%   Compute defocus strength (diopters) from distance behind lens (m) and
%   optical power of lens (diopters)
epsilon = (P.^2.*(dz./1000))./(1-P.*(dz./1000));
end

function dz = dz_from_eps(eps, P)
%   Compute dz (meters) behind lens of power P with eps diopters
%   of defocus.
f = (1/P); % Focal length (meters)
dz = (eps.*(f.^2))./((eps.*f)+1);
end

%% Local functions
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


