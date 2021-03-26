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

{% include base_path.html %}

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

The following is a list of the main families of tractography methods in
chronological order:

* Local tractography (Conturo et al. 1999, Mori et al. 1999, Basser et al.
2000).
* Global tracking (Mangin et al. 2002)
* Particle Filtering Tractography (PFT) (Girard et al. 2014)
* Parallel Transport Tractography (PTT) (Aydogan et al., 2019)

Local tractography methods and PFT can use two approaches to propagate the
streamlines:
- Deterministic: propagates streamlines consistently using the same propagation
direction.
- Probabilistic: uses a distribution function to sample from in order to decide
on the next propagation direction at each step.

Several algorithms exist to perform local tracking, depending on the
local orientation construct used or the order of the integration being
performed, among others: FACT (Mori et al. 1999), EuDX (Garyfallidis 2012),
iFOD1 (Tournier et al. 2012) / iFOD2 (Tournier et al. 2010), and SD_STREAM
(Tournier et al. 2012) are some of those. Different strategies to reduce the
uncertainty (or missed configurations) on the tracking results have also been
proposed (e.g. Ensemble Tractography (Takemura et al. 2016), Bootstrap
Tractography (Lazar et al. 2005)).

Tractography methods suffer from a number of known biases and limitations,
generally yielding tractograms containing a large number of prematurely stopped
streamlines and invalid connections, among others. This results in a hard
trade-off between sensitivity and specificity (usually measured in the form of
bundle overlap and overreach) (Maier-Hein et al. 2017).

Several enhancements to the above frameworks have been proposed, usually based
on incorporating some *a priori* knowledge (e.g. Anatomically-Constrained
Tractography (ACT) (Smith et al. 2012), Structure Tensor Informed Fiber
Tractography (STIFT) (Kleinnijenhuis et al. 2012), Surface-enhanced
Tractography (SET) (St-Onge et al. 2018), Bundle-Specific Tractography (BST)
(Rheault et al., 2019), etc.).

In the recent years, many deep learning methods have been proposed to map the
local orientation reconstruction (or directly the diffusion MRI data) to long
range white matter connectivity.
