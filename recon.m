%% recon
% Performs reconstruction using Tikhonov deconvolution for grayscale or color images.
%
% Usage:
%   recon = recon(im, psf, lambda);
%
% Inputs:
%   im      - Input image (2D grayscale or 3D color image in GRB order)
%   psf     - Point spread function for deconvolution
%   lambda  - Regularization strength(s)
%             - scalar:     Applied to all color channels 
%             - 3x1 vector: Separate values for each channel (R,G,B)
%
% Output:
%   recon   - Reconstructed image (grayscale or color RGB, same size as input)
%
% Dependencies:
%   Requires `tikhonov_deconv` for Tikhonov deconvolution.
%
% Author: Corey Simmerer
% Date  : January 6, 2025

function [recon] = recon(im, psf, lambda)

% Grayscale image
if ismatrix(im)
    recon = tikhonov_deconv(im,psf,lambda);
end

% Color image
if ~ismatrix(im)
    if ndims(im)>3
        im = squeeze(im(:,:,1,:));
    end
    
    % If only 1 regularization given, use for all 3 variables
    if length(lambda)==1
        lambda = repmat(lambda,1,3);
    end

    R = tikhonov_deconv(im(:,:,2),psf,lambda(1));
    G =  tikhonov_deconv(im(:,:,1),psf,lambda(2));
    B =  tikhonov_deconv(im(:,:,3),psf,lambda(3));

    rgb = zeros(size(im,1),size(im,2),size(im,3));
    rgb(:,:,1) = R;
    rgb(:,:,2) = G;
    rgb(:,:,3) = B;
    recon = rgb;
end

end