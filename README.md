# Computational Fundus Imaging 

This repository provides MATLAB code for computational fundus imaging and wavefront analysis using a novel diffuser-based imaging system. The approach leverages phase mask-based computational imaging to enable digital refocusing across a large range of defocus errors without requiring moving parts, as detailed in:

> **"In vivo fundus imaging and computational refocusing with a diffuser-based fundus camera"**  
> Corey Simmerer, Marisa Morakis, Lei Tian, Lia Gomez-Perez, T.Y. Alvin Liu, and Nicholas J. Durr.

## Study Highlights

- **Digital Refocusing**: Enables refocusing of fundus images computationally, covering a range of ≥ ±10 diopters of defocus error.
- **Resolution and Field of View**: Demonstrates resolutions of 7.7–9.6 line pairs per mm with a field of view ≥35 degrees.
- **In Vivo Demonstration**: Provides the first in vivo fundus images using a diffuser-based camera.
- **Ocular Safety**: Validated for human use under ISO 15004-2:2007 standards for ophthalmic devices.

This repository includes code to replicate key components of the study, such as focal stack reconstruction, resolution analysis, and optical system modeling.

## Requirements

- MATLAB R2022a or later
- Image Processing Toolbox
- Signal Processing Toolbox
- Dependencies:
  - [radialavg.zip](https://www.mathworks.com/matlabcentral/fileexchange/46468-radialavg-zip)
  - [shadedErrorBar](https://github.com/raacampbell/shadedErrorBar)

## Directory structure

```
repository/
├── tikhonov_deconv.m             # Implements Tikhonov deconvolution
├── recon.m                       # Reconstructs grayscale or color images using Tikhonov deconvolution
├── focus_stack.m                 # Refocuses raw diffuser images into a focal stack
├── resolution_analysis.m         # Analyzes system resolution by modeling with ray matrices
├── RayMat.m                      # Custom class for simulating systems with ray matrices
├── data/                         # Contains calibration files and metadata
│   ├── PSFs/                     # Contains PSF image dataset (download separately)
│   ├── RawImages/                # Contains raw diffuser images found in paper
│   ├── Output/                   # Empty, output directory for reconstructions
│   ├── calibrationparams.mat     # Calibration parameters for the imaging system
│   └── PSFinfo.mat               # Metadata for PSFs (the PSFs must be downloaded separately)
└── README.md                     # This file
```

## Components

### `tikhonov_deconv.m`

This is an implementation of Tikhonov deconvolution using the FFT.

### `recon.m`

The `recon` function reconstructs raw diffuser images using Tikhonov deconvolution and `tikhonov_deconv.m`. It supports both grayscale and color images. 

### `focus_stack.m`

This script generates a focal stack from raw diffuser images using a set of PSFs, saving it to the specified output directory. The PSFs used in this study can be downloaded from [LINK] and placed in the directory structure in `data/PSFs/`.

### `resolution_analysis.m`
Script to analyze the resolution of the system using the MTFs of the point spread function across defocus errors. The maximum theoretical resolution is defined here to be the highest spatial frequency in the radially averaged MTF that is above the noise floor of the camera. This uses the `RayMat` class to simulate the system.


### `RayMat` class

The `RayMat` class is used to simulate optical systems using ray matrices (ABCD matrices).

#### Properties

- `elements` is a cell array containing the list of optical elements in first to last order.
- `rayMatrix` is a single 2x2 ray matrix representing all elements in `elements`. It is updated upon creation of the object and when elements are added to the object.

#### Functions

- `obj = RayMat(elements)`
  - Constructor of the `RayMat` class.
  - Input: `elements`, a cell array of 2x2 matrices in first to last order.
  - Output: `RayMat` object
- `obj = addElements(elements)`
  - Appends elements to the `RayMat` system
  - Input: `elements`, a cell array of 2x2 matrices in first to last order to append
  - Output: `RayMat` object

#### Built-in optical elements

| Element                 | Function           | Parameters        |
| :----------------       | :------            | :----             |
| Free-space propagation  | `freespace(d)`     |`d`: distance      |
| Thin lens               | `thinlens(f)`      |`f`: focal length  |
| Refraction on flat surface  |  `flatRefraction(n1,n2)` | `n1`: refractive index of first medium <br/> `n2`: refractive index of second medium      |
| Refraction on spherical surface     | `curvedRefraction(n1,n2,r)`  | `n1`: refractive index of first medium <br/> `n2`: refractive index of second medium  <br> `r`: radius of curvature              |
| Thick lens    | `thickLens(n1,n2,r1,r2,t)`  | `n1`: refractive index of first medium <br/> `n2`: refractive index of second medium  <br> `r1`: radius of curvature of first surface <br> `r2`: radius of curvature of second surface   <br> `t`: Distance between surfaces           |


#### Example usage of the `RayMat` class

```matlab
% Simple 4f system
FS1 = RayMat.freespace(50);               % 50 mm propagation in free space
L1 = RayMat.thinlens(50);                 % Thin lens with focal length 50 mm
FS2 = RayMat.freespace(100);              % 100 mm propagation in free space
L2 = RayMat.thinlens(50);                 % Thin lens with focal length 50 mm
FS3 = RayMat.freespace(50);               % 50 mm propagation in free space
system = RayMat({FS1, L1, FS2, L2, FS3}); % Create the system
matrix = system.rayMatrix;                % Get the ABCD matrix of the entire system
```

## Citation

If you use this code in your work, please cite the accompanying paper:
>**"In vivo fundus imaging and computational refocusing with a diffuser-based fundus camera"**  
>Corey Simmerer, Marisa Morakis, Lei Tian, Lia Gomez-Perez, T.Y. Alvin Liu, and Nicholas J. Durr.

## Acknowledgments

This project is part of ongoing research at the Durr Lab at Johns Hopkins University, supported by The Wilmer Eye Institute PPF and the unrestricted grant to prevent blindness.

---
