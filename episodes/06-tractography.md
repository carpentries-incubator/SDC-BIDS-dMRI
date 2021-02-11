---
title: "Tractography"
teaching: 120
exercises: 20
questions:
- "What information can dMRI provide at the long range level?"
objectives:
- "Present different long range orientation reconstruction methods"
keypoints:
- "Provides an estimation of the long range underlying fiber arrangement"
- "Tractography is central to estimate and provide measures of the white matter neuroanatomy"
start: true
---

## Tractography

The local fiber orientation reconstruction can be used to map the voxel-wise
fiber orientations to white matter long range structural connectivity.
Tractography is a fiber tracking technique that studies how the local
orientations can be integrated to provide an estimation of the white matter
fibers connecting structurally two regions in the white matter.

Tractography models axonal trajectories as geometrical entities called
*streamlines* from local directional information. Tractograhy essentially uses
an integral equation involving a set of discrete local directions to numerically
find the curve (i.e. the streamline) that joins them. The streamlines generated
by a tractography method and the required meta-data are usually saved into files
called *tractograms*.

The following is a list of the main families of tractography methods:

* Local tractography (Conturo et al. 1999, Mori et al. 1999, Basser et al. 2000)
* Particle Filtering Tractography (PFT) (Girard et al. 2014)
* Global tracking (Mangin et al. 2002)

The first two methods can use two approaches to propagate the streamlines:
- Deterministic: propagates streamlines consistently using the same propagation
direction.
- Probabilistic: uses a distribution function to sample from in order to decide
on the next propagation direction at each step.

Tractography methods suffer from a number of known biases and limitations,
generally yielding tractograms containing a large number of prematurely stopped
streamlines and invalid connections, among others. This results in a hard
trade-off between sensitivity and specificity (usually measured in the form of
bundle overlap and overreach) (Maier-Hein et al. 2017).

Several enhancements to the above frameworks have been proposed (e.g. Ensemble
Tractography (Takemura et al. 2016), Surface-enhanced Tractography (SET)
(St-Onge et al. 2018), Bundle-Specific Tractography (BST) (Rheault et al.,
2019), etc.).

In the recent years, many deep learning methods have been proposed to map the
local orientation reconstruction (or directly the diffusion MRI data) to long
range white matter connectivity.
