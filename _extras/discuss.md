---
title: Discussion
---

## Setting up the camera viewpoint

The [Constrained Spherical Deconvolution] episode defines an utility function
(`generate_anatomical_views`) that receives a list of data and generates a
[Matplotlib] `Figure` instance containing the axial superior, sagittal right
and coronal anterior views of the data at the specified slices. The method
achieves this by setting up the camera in the appropriate viewpoints. These
viewpoints are obtained applying rotations at given angles; in the mentioned
method, the rotations are applied to the focal point. These rotations are then
defined according to the following perpendicular axes:

- **pitch** (transverse axis): an axis running from right to left, allowing
rotations in the sagittal plane.
- **roll** (longitudinal axis): an axis directed forward, allowing rotations in
the coronal plane.
- **yaw** (vertical axis): an axis running bottom to top, allowing rotations in
the axial plane.

The mentioned concepts are common across different 3D visualization packages.
The mentioned episode uses the [FURY] package to generate the visualization
scenes.


[Constrained Spherical Deconvolution]: https://carpentries-incubator.github.io/SDC-BIDS-dMRI/constrained_spherical_deconvolution/index.html

{% include links.md %}
