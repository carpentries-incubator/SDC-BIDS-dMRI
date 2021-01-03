---
layout: page
title: "Instructor Notes"
---
{% include base_path.html %}

## Instructor notes

## Lesson motivation and learning objectives

This lesson is designed to introduce learners to the analysis of diffusion
Magnetic Resonance Imaging (dMRI) data using primarily Python. Although it is
designed for learners who have no prior experience with diffusion MRI, a basic
understanding of MRI and basic command-line and Python skills are required.

Diffusion MRI is a feature-rich modality. Processing dMRI data involves
numerous aspects, and can take a non-negligible amount of time. Being able
to critically assess the obtained features requires practice, usually with
heterogeneous data.

Upon completion of the lesson, learners should be able to know the purpose and
value of acquiring and analyzing diffusion data, as well as being able to
process and analyze their own diffusion data. This lesson does not cover the
details about the different dMRI data acquisition sequences. Similarly, this
lesson does not cover quality assurance or checking aspects of the processing,
neither does it focus on the importance of visualization for that purpose.

## Lesson design

#### [Introduction to Diffusion MRI data]({{ relative_root_path }}/{% link _episodes/01-introduction_diffusion_data.md %})

* If your workshop includes the [Introduction to MRI and BIDS](https://carpentries-incubator.github.io/SDC-BIDS-IntroMRI/) lesson,
learners will have the necessary knowledge to better understand how diffusion
images are generated and different from other MRI data acquisitions. Similarly,
they will be aware of the reasons that prompted the neuroimaging community to
agree on a common storage convention for the data.
* Be sure that learners understand the physical phenomenon diffusion MRI is
sensitive to, and how diffusion is able to capture it.
* Explain that there are other software and tools, such as [`3D Slicer`](https://www.slicer.org/),
[`DSI Studio`](http://dsi-studio.labsolver.org/), [`ExploreDTI`](https://www.exploredti.com/),
[`MRtrix`](https://www.mrtrix.org/), or [`TrackVis`](https://www.mrtrix.org/),
to analyze or visualize diffusion MRI data.

#### [Preprocessing dMRI data]({{ relative_root_path }}/{% link _episodes/02-diffusion_preprocessing.md %})

* Pre-processing in dMRI depends on the available data (acquisition) and the
quality of the data, so learners should be encouraged to look at their data to
identify artefacts.
* Similarly, learners should embrace the importance of exploring the output at
each step in order to ensure that the quality of the data has improved after
the pre-processing step.
* Some pre-processing steps in this lesson make use of tools (such as `FSL` or
`ANTs`) that are used as command-line tools, so learners should be encouraged
to check their documentation, and adjust the arguments as necessary.

#### [Diffusion Tensor Imaging (DTI)]({{ relative_root_path }}/{% link _episodes/03-diffusion_tensor_imaging.md %})

* Learners should be able to understand the use, relevance and limitations of
the DTI model, both from the a clinical point of view, and a research setting.
* Be sure to emphasize that even if the DTI model is unable to solve fiber
crossings, the concepts and derived scalar maps are still relevant in
practice.
* Pay special attention to the visualization of the DTI model results in order
to ensure that the anatomical and the diffusion MR images are correctly
registered.
* Also, stress the fact that DTI is a model applied to diffusion MRI, and that
DTI must not be mistaken with the modality itself (with the acronyms DWI and
DTI mistakenly being used interchangeably, partially due to the extensive use
of the DTI model in clinical practice).

#### [Tractography]({{ relative_root_path }}/{% link _episodes/04-tractography.md %})

* Make sure to explain the difference between the actual biological white matter
fibers and streamlines in tractograms, and why tractography is not quantitative
in terms of connectivity.
* Emphasize the limitations of tractography in downstream tasks, and the need
to post-process tractograms to reliably account fo the white matter anatomy.
* It may be necessary to let the learners know about the potential
inconsistencies in tractography file formats. These may be due to different
conventions when choosing the origin or center across tools. At times it is
necessary to visualize the tractograms using different software.

#### Concluding remarks

* Try to provide the learners with a clear idea of where diffusion MRI sits
among the rest of MRI or neuroimaging modalities.
* Make sure that all learners are able to follow the lessons and acomodate
the pace to the audience. Diffusion MRI being such a feature-rich imaging
modality, it is important that the explained concepts are clear to the
learners before continuing to the next aspect. Be prepared to be flexible and
skip some aspects if necessary.
* Although the exercises are concise, be generous in the time allocated so
that if part of the audience requires additional explanations, helpers can
assist them during that time.
* When possible, provide remarks about avenues that diffusion MRI researchers
may be exploring to gain further insight from dMRI data and overcome the
limitations of current methods.

## Technical tips and tricks

* Be clear about the purpose and convenience of using `Jupyter Notebooks` to
teach the lesson. Allow some time at the beginning of the day to provide an
overview of how using the learned tools would be translated to a more formal
dMRI data analysis setting once a research aspect or analysis has been
consolidated.
* Provide a broad overview of the software tools and packages used throughout
the lesson, and the convenience of each of them.
* Be sure to link the tools with their corresponding documentation pages.
Being at ease reading the documentation and understanding what parameters
mean will be essential when learners will analyze their own data back at
their institution or research works.
* If time permits, use mistakes as teaching moments: the most vital skill
you can impart is how to debug and recover from unexpected errors.

## Common problems




{% include links.md %}
