---
title: "Introduction to Diffusion"
teaching: 20
exercises: 5
questions:
- "How is dMRI data represented stored"
- "What is diffusion weighting"
objectives:
- "Representation of diffusion data and associated gradients"
- "Learn about diffusion gradients"
keypoints:
- "dMRI data is represented as a 4-dimensional image (x,y,z,diffusion directional sensitivity)"
- "dMRI data is sensitive to a particular direction of diffusion motion. Due to this sensitivity, each volume of the 4D image is sensitive to a particular direction"
---

## Diffusion Weighted Imaging (DWI)

Diffusion imaging probes the random, microscopic motion of water protons by employing MRI sequences which are sensitive to the geometry and environmental organization surrounding the water protons. This is a popular technique for studying the white matter of the brain. The diffusion within biological structures, such as the brain, are often restricted due to barriers (eg. cell membranes), resulting in a preferred direction of diffusion (anisotropy). A typical diffusion MRI scan will acquire multiple volumes that are sensitive to a particular diffusion direction and result in diffusion-weighted images (DWI). Diffusion that exhibits directionality in the same direction result in an attenuated signal. With further processing (to be discussed later in the lesson), the acquired images can provide measurements which are related to the microscopic changes and estimate white matter trajectories. Images with no diffusion weighting are also acquired as part of the acquisition protocol.

![fiber_configurations](../fig/DiffusionDirections.png){:class="img-responsive"} \
Diffusion along X, Y, and Z directions

## b-values & b-vectors

In addition to the acquired diffusion images, two files are collected as part of the diffusion dataset. These files correspond to the gradient amplitude (b-values) and directions (b-vectors) of the diffusion measurement and are named with the extensions <code>.bval</code> and <code>.bvec</code> respectively. The b-value is the diffusion-sensitizing factor, and reflects the timing & strength of the gradients used to acquire the diffusion-weighted images. The b-vector corresponds to the direction of the diffusion sensitivity. Together these two files define the diffusion MRI measurement as a set of gradient directions and corresponding amplitudes.

## Dataset

In addition to the acquired diffusion images, two files are collected as part of the diffusion dataset. These files correspond to the gradient amplitude (b-values) and directions (b-vectors) of the diffusion measurement and are named with the extensions <code>.bval and <code>.bvec</code> respectively. The b-value is the diffusion-sensitizing factor, and reflects the timing & strength of the gradients used to acquire the diffusion-weighted images. The b-vector corresponds to the direction of the diffusion sensitivity. Together these two files define the diffusion MRI measurement as a set of gradient directions and corresponding amplitudes.

Below is a tree diagram showing the folder structure of a single MR session within ds000030. This was obtained by using the bash command <code>tree<code>.  

~~~
tree '../data/ds000030'
~~~
{: .language-bash}

~~~
../data/ds000030
├── CHANGES
├── code
├── dataset_description.json
├── README
└── sub-10788/
    ├── anat
    │   ├── sub-10788_T1w.json
    │   └── sub-10788_T1w.nii.gz
    └── dwi
        ├── sub-10788_dwi.bval
        │── sub-10788_dwi.bvec
        │── sub-10788_dwi.json
        └── sub-10788_dwi.nii.gz
~~~
{: .output}

## Querying a BIDS Dataset

[`pybids`](https://bids-standard.github.io/pybids/) is a Python API for querying, summarizing and manipulating the BIDS folder structure. We will make use of <code>pybids</code> to query the necessary files.

Lets first pull the metadata from its associated JSON file using the <code>get_metadata()</code> function for the first run.

~~~
from bids.layout import BIDSLayout

layout = BIDSLayout("../../data/ds000030", validate=False)
~~~
{: .language-python}

Now that we have a layout object, we can work with a BIDS dataset! Lets extract the metadata. from the dataset.

~~~
dwi = layout.get(subject='10788', suffix='dwi', extension='nii.gz', return_type='file')[0]
layout.get_metadata(dwi)
~~~
{: .language-python}

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
dwi = layout.get(subject='10788', suffix='dwi', extension='.nii.gz', return_type='file')[0]
bvec = layout.get(subject='10788', suffix='dwi', extension='bvec', return_type='file')[0]
bval = layout.get(subject='10788', suffix='dwi', extension='bval', return_type='file')[0]
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
(96, 96, 60, 65)
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

![DWI Slice](../fig/dwi_slice.png){:class="img-responsive"}

We can also see how the diffusion gradients are represented. This is plotted on a sphere, the further away from the center of the sphere, the stronger the diffusion gradient (increased sensitivity to diffusion).

~~~
bvec_txt = np.genfromtxt(bvec)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(bvec_txt[0], bvec_txt[1], bvec_txt[2])
~~~
{: .language-python}

![Diffusion Gradient Sphere](../fig/diffusion_gradient.png){:class="img-responsive"}

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
array([[-9.99984e-01,  4.03613e-03,  4.03613e-03],
       [-9.92980e-04, -9.99987e-01,  4.98886e-03],
       [-2.48897e-02,  6.53243e-01, -7.56739e-01],
       [-5.89518e-01,  7.69236e-01,  2.46462e-01],
       [ 2.35626e-01,  5.28739e-01,  8.15423e-01],
       [ 8.93067e-01,  2.63666e-01,  3.64570e-01],
       [ 7.97398e-01,  1.33552e-01, -5.88489e-01],
       [ 2.32919e-01,  9.31812e-01, -2.78344e-01],
       [ 9.36380e-01,  1.44036e-01, -3.20072e-01],
       [-5.04032e-01,  8.46814e-01, -1.69873e-01],
       [-3.44841e-01,  8.50410e-01, -3.97351e-01],
       [-4.55850e-01,  6.35469e-01, -6.23202e-01],
       [ 4.87386e-01,  3.93024e-01,  7.79735e-01],
       [-6.16792e-01,  6.76545e-01, -4.02310e-01],
       [ 5.77851e-01,  1.09487e-01, -8.08765e-01],
       [ 8.25555e-01,  5.24662e-01,  2.07818e-01],
       [ 8.94898e-01, -4.48150e-02,  4.44016e-01],
       [-2.89332e-01,  5.45724e-01, -7.86430e-01],
       [-1.15014e-01,  9.64050e-01, -2.39541e-01],
       [ 8.00058e-01, -4.08056e-01, -4.39770e-01],
       [ 5.11970e-01,  8.42290e-01, -1.68625e-01],
       [ 7.89764e-01, -1.57178e-01, -5.92932e-01],
       [ 9.49115e-01, -2.37601e-01, -2.06705e-01],
       [ 2.32032e-01,  7.86655e-01, -5.72132e-01],
       [ 1.96515e-02,  1.91844e-01, -9.81229e-01],
       [-2.15961e-01,  9.57087e-01,  1.93247e-01],
       [ 7.72435e-01, -6.07408e-01, -1.85471e-01],
       [-1.59879e-01,  3.59797e-01, -9.19231e-01],
       [-1.46103e-01,  7.34950e-01,  6.62195e-01],
       [ 8.87180e-01,  4.21444e-01, -1.87872e-01],
       [-5.62338e-01,  2.36544e-01,  7.92353e-01],
       [-3.80669e-01,  1.46788e-01, -9.12987e-01],
       [ 3.05803e-01,  2.02751e-01, -9.30256e-01],
       [ 3.32094e-01,  1.33876e-01,  9.33697e-01],
       [ 9.62206e-01,  2.69443e-01,  3.95029e-02],
       [ 9.59295e-01, -2.09888e-01,  1.88943e-01],
       [-4.50964e-01,  8.90337e-01,  6.27015e-02],
       [ 7.71192e-01, -6.31175e-01,  8.29533e-02],
       [-7.09223e-01, -4.12894e-01,  5.71421e-01],
       [-6.94205e-01,  2.78961e-02, -7.19236e-01],
       [ 6.81181e-01,  5.33350e-01,  5.01528e-01],
       [-1.40978e-01, -7.29050e-01, -6.69784e-01],
       [ 7.40351e-01, -3.93222e-01,  5.45212e-01],
       [-1.01944e-01,  8.25404e-01, -5.55261e-01],
       [ 5.83509e-01, -6.00385e-01, -5.46859e-01],
       [ 8.66669e-02,  3.39104e-01,  9.36748e-01],
       [ 5.50506e-01,  7.95484e-01,  2.53276e-01],
       [ 8.37371e-01, -4.62163e-01,  2.91916e-01],
       [-3.62527e-01,  5.65304e-01,  7.40949e-01],
       [-1.83461e-01,  3.96756e-01,  8.99404e-01],
       [-7.18319e-01, -6.95701e-01, -4.24514e-03],
       [ 4.31996e-01,  6.86464e-01,  5.84933e-01],
       [ 5.00977e-01,  6.94308e-01, -5.16680e-01],
       [ 1.69597e-01,  5.13550e-01, -8.41132e-01],
       [ 4.63360e-01,  4.27481e-01, -7.76246e-01],
       [ 3.84024e-01, -8.12297e-01, -4.38975e-01],
       [ 7.13857e-01,  2.51359e-01,  6.53626e-01],
       [ 2.58398e-01,  8.87277e-01,  3.82061e-01],
       [-9.28270e-04,  8.02399e-02,  9.96775e-01],
       [-3.63633e-02,  9.04616e-01,  4.24675e-01],
       [-5.70681e-01,  3.07326e-01, -7.61495e-01],
       [-2.82028e-01,  1.48741e-01,  9.47806e-01],
       [ 7.19926e-01,  6.12166e-01, -3.27047e-01],
       [ 2.65067e-01,  9.60908e-01,  7.99761e-02]])
~~~
{: .output}

It is also important to know where our diffusion weighting free measurements are as we need them for registration in our preprocessing, (our next notebook). The gtab.b0s_mask shows that this is our first volume of our dataset.

~~~
gtab.b0s_mask
~~~
{: .language-python}

~~~
array([ True, False, False, False, False, False, False, False, False,
       False, False, False, False, False, False, False, False, False,
       False, False, False, False, False, False, False, False, False,
       False, False, False, False, False, False, False, False, False,
       False, False, False, False, False, False, False, False, False,
       False, False, False, False, False, False, False, False, False,
       False, False, False, False, False, False, False, False, False,
       False, False])
~~~
{: .output}

In the next few notebooks, we will talk more about preprocessing the diffusion weighted images and reconstructing the Tensor model

> ## Exercise 1
> 1. Get a list of **all** diffusion data in Nifti file format
> 2. Get the metadata for the diffusion associated with subject 10788
> > ## Solution
> > ## *All the diffusion data*
> > ~~~
> > dwi_data = layout.get('suffix=dwi', extension='nii.gz', return_type='file')
> > ~~~
> > {: .language-python}
> >
> > ## *Metadata for subject 10788*
> > ~~~
> > dwi = layout.get(subject='10788', suffix='dwi', extension='nii.gz', return_type='file')[0]
> > layout.get_metadata(dwi)
> > ~~~
> > {: .language-python}
> {: .solution}
{: .challenge}

{% include links.md %}
