---
title: "Introduction to Diffusion MRI data"
teaching: 20
exercises: 5
questions:
- "How is dMRI data represented?"
- "What is diffusion weighting?"
objectives:
- "Representation of diffusion data and associated gradients"
- "Learn about diffusion gradients"
keypoints:
- "dMRI data is represented as a 4-dimensional image (x,y,z,diffusion directional sensitivity)"
- "dMRI data is sensitive to a particular direction of diffusion motion. Due to this sensitivity, each volume of the 4D image is sensitive to a particular direction"
start: true
---

## Diffusion Weighted Imaging (DWI)

Diffusion imaging probes the random, microscopic motion of water protons by employing MRI sequences which are sensitive to the geometry and environmental organization surrounding the water protons. This is a popular technique for studying the white matter of the brain. The diffusion within biological structures, such as the brain, are often restricted due to barriers (eg. cell membranes), resulting in a preferred direction of diffusion (anisotropy). A typical diffusion MRI scan will acquire multiple volumes that are sensitive to a particular diffusion direction and result in diffusion-weighted images (DWI). Diffusion that exhibits directionality in the same direction result in an attenuated signal. With further processing (to be discussed later in the lesson), the acquired images can provide measurements which are related to the microscopic changes and estimate white matter trajectories. Images with no diffusion weighting are also acquired as part of the acquisition protocol.

![fiber_configurations](../fig/introduction/DiffusionDirections.png){:class="img-responsive"} \
Diffusion along X, Y, and Z directions

## b-values & b-vectors

In addition to the acquired diffusion images, two files are collected as part of the diffusion dataset. These files correspond to the gradient amplitude (b-values) and directions (b-vectors) of the diffusion measurement and are named with the extensions <code>.bval</code> and <code>.bvec</code> respectively. The b-value is the diffusion-sensitizing factor, and reflects the timing & strength of the gradients used to acquire the diffusion-weighted images. The b-vector corresponds to the direction of the diffusion sensitivity. Together these two files define the diffusion MRI measurement as a set of gradient directions and corresponding amplitudes.

## Dataset

In addition to the acquired diffusion images, two files are collected as part of the diffusion dataset. These files correspond to the gradient amplitude (b-values) and directions (b-vectors) of the diffusion measurement and are named with the extensions <code>.bval</code> and <code>.bvec</code> respectively. The b-value is the diffusion-sensitizing factor, and reflects the timing & strength of the gradients used to acquire the diffusion-weighted images. The b-vector corresponds to the direction of the diffusion sensitivity. Together these two files define the diffusion MRI measurement as a set of gradient directions and corresponding amplitudes.

Below is a tree diagram showing the folder structure of a single MR subject and session within ds000221. This was obtained by using the bash command <code>tree<code>.

~~~
tree '../data/ds000221'
~~~
{: .language-bash}

~~~
../data/ds000221
├── .bidsignore
├── CHANGES
├── dataset_description.json
├── participants.tsv
├── README
├── derivatives/
├── sub-010001/
└── sub-010002/
    ├── ses-01/
    │    ├── anat
    │    │    ├── sub-010002_ses-01_acq-lowres_FLAIR.json
    │    │    ├── sub-010002_ses-01_acq-lowres_FLAIR.nii.gz
    │    │    ├── sub-010002_ses-01_acq-mp2rage_defacemask.nii.gz
    │    │    ├── sub-010002_ses-01_acq-mp2rage_T1map.nii.gz
    │    │    ├── sub-010002_ses-01_acq-mp2rage_T1w.nii.gz
    │    │    ├── sub-010002_ses-01_inv-1_mp2rage.json
    │    │    ├── sub-010002_ses-01_inv-1_mp2rage.nii.gz
    │    │    ├── sub-010002_ses-01_inv-2_mp2rage.json
    │    │    ├── sub-010002_ses-01_inv-2_mp2rage.nii.gz
    │    │    ├── sub-010002_ses-01_T2w.json
    │    │    └── sub-010002_ses-01_T2w.nii.gz
    │    ├── dwi
    │    │    ├── sub-010002_ses-01_dwi.bval
    │    │    │── sub-010002_ses-01_dwi.bvec
    │    │    │── sub-010002_ses-01_dwi.json
    │    │    └── sub-010002_ses-01_dwi.nii.gz
    │    ├── fmap
    │    │    ├── sub-010002_ses-01_acq-GEfmap_run-01_magnitude1.json
    │    │    ├── sub-010002_ses-01_acq-GEfmap_run-01_magnitude1.nii.gz
    │    │    ├── sub-010002_ses-01_acq-GEfmap_run-01_magnitude2.json
    │    │    ├── sub-010002_ses-01_acq-GEfmap_run-01_magnitude2.nii.gz
    │    │    ├── sub-010002_ses-01_acq-GEfmap_run-01_phasediff.json
    │    │    ├── sub-010002_ses-01_acq-GEfmap_run-01_phasediff.nii.gz
    │    │    ├── sub-010002_ses-01_acq-SEfmapBOLDpost_dir-AP_epi.json
    │    │    ├── sub-010002_ses-01_acq-SEfmapBOLDpost_dir-AP_epi.nii.gz
    │    │    ├── sub-010002_ses-01_acq-SEfmapBOLDpost_dir-PA_epi.json
    │    │    ├── sub-010002_ses-01_acq-SEfmapBOLDpost_dir-PA_epi.nii.gz
    │    │    ├── sub-010002_ses-01_acq-sefmapBOLDpre_dir-AP_epi.json
    │    │    ├── sub-010002_ses-01_acq-sefmapBOLDpre_dir-AP_epi.nii.gz
    │    │    ├── sub-010002_ses-01_acq-sefmapBOLDpre_dir-PA_epi.json
    │    │    ├── sub-010002_ses-01_acq-sefmapBOLDpre_dir-PA_epi.nii.gz
    │    │    ├── sub-010002_ses-01_acq-SEfmapBOLDpost_dir-AP_epi.json
    │    │    ├── sub-010002_ses-01_acq-SEfmapBOLDpost_dir-AP_epi.nii.gz
    │    │    ├── sub-010002_ses-01_acq-SEfmapBOLDpost_dir-PA_epi.json
    │    │    ├── sub-010002_ses-01_acq-SEfmapBOLDpost_dir-PA_epi.nii.gz
    │    │    ├── sub-010002_ses-01_acq-SEfmapDWI_dir-AP_epi.json
    │    │    ├── sub-010002_ses-01_acq-SEfmapDWI_dir-AP_epi.nii.gz
    │    │    ├── sub-010002_ses-01_acq-SEfmapDWI_dir-PA_epi.json
    │    │    └── sub-010002_ses-01_acq-SEfmapDWI_dir-PA_epi.nii.gz
    │    └── fmap
    │    │    ├── sub-010002_ses-01_task-rest_acq-AP_run-01_bold.json
    │    │    └── sub-010002_ses-01_task-rest_acq-AP_run-01_bold.nii.gz
    └── ses-02/
~~~
{: .output}

## Querying a BIDS Dataset

[`pybids`](https://bids-standard.github.io/pybids/) is a Python API for querying, summarizing and manipulating the BIDS folder structure. We will make use of <code>pybids</code> to query the necessary files.

Lets first pull the metadata from its associated JSON file using the <code>get_metadata()</code> function for the first run.

~~~
from bids.layout import BIDSLayout

layout = BIDSLayout("../../data/ds000221", validate=False)
~~~
{: .language-python}

Now that we have a layout object, we can work with a BIDS dataset! Lets extract the metadata. from the dataset.

~~~
dwi = layout.get(subject='010002', suffix='dwi', extension='nii.gz', return_type='file')[0]
layout.get_metadata(dwi)
~~~
{: .language-python}

~~~
{'EchoTime': 0.08,
 'EffectiveEchoSpacing': 0.000390001,
 'FlipAngle': 90,
 'ImageType': ['ORIGINAL', 'PRIMARY', 'DIFFUSION', 'NON'],
 'MagneticFieldStrength': 3,
 'Manufacturer': 'Siemens',
 'ManufacturersModelName': 'Verio',
 'MultibandAccelerationFactor': 2,
 'ParallelAcquisitionTechnique': 'GRAPPA',
 'ParallelReductionFactorInPlane': 2,
 'PartialFourier': '7/8',
 'PhaseEncodingDirection': 'j-',
 'RepetitionTime': 7,
 'TotalReadoutTime': 0.04914}
~~~
{: .output}

## Diffusion Imaging in Python ([dipy](https://dipy.org))

For this lesson, we will use the <code>Dipy</code> (Diffusion Imaging in Python) package for processing and analysing diffusion MRI.

### Why <code>dipy</code>?

- Fully free and open source
- Implemented in Python. Easy to understand, and easy to use.
- Implementations of many state-of-the art algorithms
- High performance. Many algorithms implemented in [cython](http://cython.org/)

### Installing <code>dipy</code>

The easiest way to install <code>Dipy</code> is to use <code>pip</code>! Additionally, <code>Dipy</code> makes use of the FURY library for visualization. We will also install this using <code>pip</code>!

We can install it by entering the following in a terminal <code>pip install dipy</code>. We will do so using Jupyter Magic in the following cell!

### Defining a measurement: <code>GradientTable</code>

<code>Dipy</code> has a built-in function that allows us to read in <code>bval</code> and <code>bvec</code> files named <code>read_bvals_bvecs</code> under the <code>dipy.io.gradients</code> module. Let's first grab the path to our gradient directions and amplitude files and load them into memory.

~~~
dwi = layout.get(subject='010002', suffix='dwi', extension='.nii.gz', return_type='file')[0]
bvec = layout.get(subject='010002', suffix='dwi', extension='bvec', return_type='file')[0]
bval = layout.get(subject='010002', suffix='dwi', extension='bval', return_type='file')[0]
~~~
{: .language-python}

Now that we have the necessary diffusion files, lets explore the data!

~~~
import numpy as np
import nibabel as nib

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

data = nib.load(dwi).get_fdata()
data.shape
~~~
{: .language-python}

~~~
(128, 128, 88, 67)
~~~
{: .output}

We can see that the data is 4 dimensional. The 4th dimension represents the different diffusion directions we are sensitive to. Next, let's take a look at a slice.

~~~
x_slice = data[58, :, :, 0]
y_slice = data[:, 58, :, 0]
z_slice = data[:, :, 30, 0]

slices = [x_slice, y_slice, z_slice]

fig, axes = plt.subplots(1, len(slices))
for i, slice in enumerate(slices):
    axes[i].imshow(slice.T, cmap="gray", origin="lower")
~~~
{: .language-python}

![DWI Slice](../fig/introduction/dwi_slice.png){:class="img-responsive"}

We can also see how the diffusion gradients are represented. This is plotted on a sphere, the further away from the center of the sphere, the stronger the diffusion gradient (increased sensitivity to diffusion).

~~~
bvec_txt = np.genfromtxt(bvec)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(bvec_txt[0], bvec_txt[1], bvec_txt[2])
~~~
{: .language-python}

![Diffusion Gradient Sphere](../fig/introduction/diffusion_gradient.png){:class="img-responsive"}

The files associated with the diffusion gradients need to converted to a <code>GradientTable</code> object to be used with <code>Dipy</code>. A <code>GradientTable</code> object can be implemented using the <code>dipy.core.gradients</code> module. The input to the <code>GradientTable</code> should be our the values for our gradient directions and amplitudes we read in.

~~~
from dipy.io.gradients import read_bvals_bvecs
from dipy.core.gradients import gradient_table

gt_bvals, gt_bvecs = read_bvals_bvecs(bval, bvec)
gtab = gradient_table(gt_bvals, gt_bvecs)
~~~
{: .language-python}

We will need this gradient table later on to process our data and generate diffusion tensor images (DTI)!

There is also a built in function for gradient tables, <code>b0s_mask</code> that can be used to separate diffusion weighted measurements from non-diffusion weighted measurements (b=0s/mm^2). We will extract the vector corresponding to diffusion weighted measurements!

~~~
gtab.bvecs[~gtab.b0s_mask]
~~~
{: .language-python}

~~~
array([[ 0.0480948 , -0.518981  ,  0.853432  ],
       [ 0.980937  , -0.0827268 , -0.175836  ],
       [-0.24275   , -0.355443  , -0.902625  ],
       [-0.292642  , -0.878897  ,  0.376696  ],
       [ 0.085518  , -0.362038  , -0.928232  ],
       [ 0.470646  , -0.695075  ,  0.543473  ],
       [-0.865701  , -0.485398  ,  0.122274  ],
       [ 0.226775  , -0.23832   ,  0.944339  ],
       [ 0.334443  , -0.89435   , -0.29713   ],
       [ 0.727534  ,  0.0823542 ,  0.681111  ],
       [-0.625823  ,  0.126479  ,  0.769642  ],
       [ 0.353667  , -0.886595  ,  0.298111  ],
       [-0.853084  , -0.472587  , -0.221153  ],
       [ 0.516586  , -0.856143  , -0.0125875 ],
       [-0.766369  , -0.458916  ,  0.449527  ],
       [-0.125754  , -0.637058  , -0.760489  ],
       [-0.149251  , -0.743085  ,  0.652341  ],
       [ 0.937341  , -0.348404  , -0.00268246],
       [-0.124163  , -0.219376  ,  0.967708  ],
       [-0.36458   , -0.906103  , -0.214613  ],
       [ 0.579767  , -0.204866  ,  0.788606  ],
       [ 0.445586  , -0.692181  , -0.567749  ],
       [-0.905294  , -0.174158  , -0.387442  ],
       [-0.282606  , -0.482101  ,  0.829284  ],
       [-0.731609  , -0.431568  , -0.527729  ],
       [ 0.676757  , -0.438047  , -0.591705  ],
       [ 0.171374  , -0.758177  ,  0.629125  ],
       [-0.59121   , -0.711294  , -0.380173  ],
       [-0.46017   , -0.155292  ,  0.874144  ],
       [ 0.845645  , -0.162694  ,  0.508345  ],
       [-0.130717  , -0.98886   ,  0.0712005 ],
       [ 0.975624  , -0.0906642 ,  0.199845  ],
       [-0.288147  ,  0.100936  ,  0.952252  ],
       [ 0.655193  , -0.70285   ,  0.276989  ],
       [-0.442479  , -0.607006  , -0.660118  ],
       [-0.471845  , -0.674107  ,  0.56828   ],
       [ 0.638596  , -0.702848  , -0.313369  ],
       [-0.642432  , -0.712214  ,  0.2829    ],
       [ 0.850936  , -0.431779  ,  0.299125  ],
       [-0.240808  , -0.830023  , -0.503063  ],
       [-0.578162  , -0.401657  ,  0.710211  ],
       [ 0.100487  , -0.838741  , -0.535178  ],
       [-0.924592  , -0.196278  ,  0.326504  ],
       [-0.0210952 , -0.967749  , -0.251032  ],
       [-0.764669  , -0.142405  ,  0.628492  ],
       [ 0.197294  , -0.980191  ,  0.0173175 ],
       [ 0.405727  ,  0.0420623 ,  0.913026  ],
       [ 0.859032  , -0.144763  , -0.491028  ],
       [ 0.380277  , -0.486027  ,  0.786872  ],
       [-0.6891    , -0.721525  , -0.0674049 ],
       [ 0.430722  , -0.388461  , -0.814602  ],
       [ 0.0366712 , -0.92944   ,  0.367147  ],
       [-0.540564  , -0.318621  , -0.778634  ],
       [ 0.775224  , -0.631369  , -0.0200052 ],
       [ 0.0646129 ,  0.0600214 ,  0.996104  ],
       [-0.978577  , -0.203636  , -0.0303294 ],
       [ 0.199971  , -0.618334  , -0.760049  ],
       [ 0.678143  , -0.446978  ,  0.583381  ],
       [-0.448761  , -0.888954  ,  0.0915084 ],
       [ 0.849148  , -0.426713  , -0.311228  ]])
~~~
{: .output}

It is also important to know where our diffusion weighting free measurements are as we need them for registration in our preprocessing, (our next notebook). The gtab.b0s_mask shows that this is our first volume of our dataset.

~~~
gtab.b0s_mask
~~~
{: .language-python}

~~~
array([ True, False, False, False, False, False, False, False, False,
       False, False,  True, False, False, False, False, False, False,
       False, False, False, False,  True, False, False, False, False,
       False, False, False, False, False, False,  True, False, False,
       False, False, False, False, False, False, False, False,  True,
       False, False, False, False, False, False, False, False, False,
       False,  True, False, False, False, False, False, False, False,
       False, False, False,  True])
~~~
{: .output}

In the next few notebooks, we will talk more about preprocessing the diffusion weighted images, reconstructing the diffusion tensor model, and reconstruction axonal trajectories via tractography.

> ## Exercise 1
> 1. Get a list of **all** diffusion data in Nifti file format
> 2. Get the metadata for the diffusion associated with subject 010016.
> > ## Solution
> > ## *All the diffusion data for subject 010016*
> > ~~~
> > dwi_data = layout.get(suffix=dwi', extension='nii.gz', return_type='file')
> > ~~~
> > {: .language-python}
> >
> > ## *Metadata for subject 010016*
> > ~~~
> > dwi = layout.get(subject='010016', suffix='dwi', extension='nii.gz', return_type='file')[0]
> > layout.get_metadata(dwi)
> > ~~~
> > {: .language-python}
> {: .solution}
{: .challenge}

{% include links.md %}
