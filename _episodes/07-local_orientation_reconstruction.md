---
title: "Local fiber orientation reconstruction"
teaching: 120
exercises: 20
questions:
- "What information can dMRI provide at the voxel level?"
objectives:
- "Present different local orientation reconstruction methods"
keypoints:
- "Provides an estimation of the local (voxel-wise) underlying fiber orientation"
- "Local fiber orientation reconstruction is the primer to all dMRI derivatives"
---

## Orientation reconstruction

Diffusion MRI is sensitive to the underlying white matter fiber orientation
distribution. Once the data has been pre-processed to remove noise and other
acquisition artefacts, dMRI data can be used to extract features that describe
the white matter. Estimating or reconstructing the local fiber orientation is
the first step to gain such insight.

The local fiber orientation reconstruction task faces several challenges
derived, among others, by the orientational heterogeneity that the white matter
presents at the voxel level. Due to the arrangement of the white matter fibers,
and the scale and limits of the diffusion modality itself, a large amount of
voxels are traversed by several fibers. Resolving such configurations with
incomplete information is not a solved task. Several additional factors, such as
imperfect models or their choices, influence the reconstruction results, and
hence the downstream results.

The following is a (non-exhaustive) list of the known local orientation
reconstruction methods:

* Diffusion Tensor Imaging (DTI) (Basser et al., 1994)
* Q-ball Imaging (QBI) (Tuch 2004; Descoteaux et al. 2007)
* Diffusion Kurtosis Imaging (DKI) (Jensen et al. 2005)
* Constrained Spherical Deconvolution (CSD) (Tournier et al. 2007)
* Diffusion Spectrum Imaging (DSI) (Wedeen et al. 2008)
* Simple Harmonic Oscillator based Reconstruction and Estimation (SHORE) (Özarslan et al. 2008)
* Constant Solid Angle (CSA) (Aganj et al. 2009)
* Damped Richardson-Lucy Spherical Deconvolution (dRL-SD) (Dell'Acqua et al. 2010)
* DSI with deconvolution (Canales-Rodriguez et al. 2010)
* Generalized Q-sampling Imaging (Yeh et al. 2010)
* Orientation Probability Density Transform (OPDT) (Tristan-Vega et al. 2010)
* Mean Apparent Propagator (MAPMRI) (Özarslan et al. 2013)
* Sparse Fascicle Model (SFM) (Rokem et al. 2015)
* Robust and Unbiased Model-Based Spherical Deconvolution (RUMBA-SD) (Canales-Rodriguez et al. 2015)
* Sparse Bayesian Learning (SBL) (Canales-Rodriguez et al. 2019)

These methods vary in terms of the required data. Hence, there are a few factors
that influence the choice for a given reconstruction method:
- The available data in terms of the number of (b-value) shells (single- or
multi-shell).
- The acquisition/sampling scheme.
- The available time to reconstruct the data.

Besides such requirements, the preference over a method generally lies in its
ability to resolve complex fiber configurations, such as fiber crossings at
reduced angles. Additionally, some of these methods provide additional products
beyond the orientation reconstruction that might also be of interest.

Finally, several deep learning-based methods have been proposed to estimate the
local fiber orientation using the diffusion MRI data.
