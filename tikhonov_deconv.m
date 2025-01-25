%% tikhonov_deconv
% This function implements Tikhonov deconvolution via the regularized
% inverse filter using FFTs.
%
% Usage: 
%   x_hat = tikhonov_deconv(y, h, lambda, removeBackground=True);
%
% Inputs:
%   y       - Raw input image (2D array normalized between 0 and 1)
%   h       - Point spread function (2D array)
%   lambda  - Regularization strength (e.g., 0.007)
%   varagin - Name-value pairs for optional parameters
%       'removeBackground' - Subtracts the average value from top right
%       corner of image. (Default: true)
%       'padArrays' - Zero pads images before computing FFT. (Default:
%       false)
%   
% Output:
%   x_hat   - Deconvolved image (normalized between 0 and 1)

function x_hat = tikhonov_deconv(y, h, lambda, varargin)
p = inputParser;
addOptional(p,'removeBackground', true, @isLogical);
addOptional(p,'padArrays', false, @isLogical);
parse(p, varargin{:});
%% 1: Remove background (optional)
y = double(y);
if p.Results.removeBackground
    y_avg = mean(y(1:50,1:50),'all'); 
    y = y - y_avg;
end
%% 2: Pad arrays (optional)
if p.Results.padArrays
    yp = padarray(y,[size(h,1),size(h,2)],"circular","both");
    hp = padarray(h, [size(y,1), size(y,2)],0, "both");
else
    yp = y;
    hp = h;
end
%% 3: Compute reconstruction
% Fourier transforms
Yp = fft2(yp); Hp = fft2(hp);

% Normalize lambda by max(abs(Hp))
r = lambda*max(abs(Hp(:)));

% Compute X_hat by regularized deconvolution
X_hat = (Yp.*conj(Hp))./(abs(Hp).^2+r.^2);

% Convert back to spatial domain and keep positive, real component
x_hat = real(ifft2(X_hat));
x_hat(x_hat<0)=0;
x_hat = fftshift(x_hat(1:size(y,1), 1:size(y,2)));

% Normalize between 0 and 1
x_hat = x_hat./max(x_hat(:));

end
