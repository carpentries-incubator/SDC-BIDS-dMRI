# Introduction to dMRI

[![Build Status](https://github.com/carpentries-incubator/SDC-BIDS-dMRI/workflows/Build,%20test/badge.svg)](https://github.com/carpentries-incubator/SDC-BIDS-dMRI/actions?query=workflow%3A"Build%2C+test")
[![Create a Slack Account with us][create_slack_svg]][slack_heroku_invite]
[![Slack Status][slack_channel_status]][slack_channel_url]
[![Binder][binder_svg]][binder_url]

An introduction to diffusion Magnetic Resonance Imaging (dMRI) analysis in
Python.

## Why Python?

Python is rapidly becoming the standard language for data analysis, visualization and automated workflow building. It is a free and open-source software that is relatively easy to pick up by new programmers. In addition, with Python packages such as `Jupyter` one can keep an interactive code journal of analysis - this is what we'll be using in the workshop. Using Jupyter notebooks allows you to keep a record of all the steps in your analysis, enabling transparency and ease of code sharing.

Another advantage of Python is that it is maintained by a large user-base. Anyone can easily make their own Python packages for others to use. Therefore, there exists a *very* large codebase for you to take advantage of for your neuroimaging analysis; from basic statistical analysis, to brain visualization tools, to advanced machine learning and multivariate methods!

## About the Lesson

This lesson teaches:
- What diffusion Magnetic Resonance Imaging is
- How dMRI data is organized within the BIDS framework
- What the standard preprocessing steps in dMRI are
- How local fiber orientation can be reconstructed using dMRI data
- How dMRI can provide insight into structural white matter connectivity

## Episodes

|   Topic  | Time | Episode | Question(s) |
|:---------|:----:|:--------|:------------|
| **Introduction to Diffusion MRI data** | 30 | [1 Introduction to Diffusion MRI data][episode01] | How is dMRI data represented?<br />What is diffusion weighting? |
| **Preprocessing dMRI data** | 30 | [2 Preprocessing dMRI data][episode02] | What are the standard preprocessing steps?<br />How do we register with an anatomical image? |
| **Local fiber orientation reconstruction** | 30 | [3 Local fiber orientation reconstruction][episode03] | What information can dMRI provide at the voxel level? |
| | 30 | [3.1 Diffusion Tensor Imaging (DTI)][episode04] | What is diffusion tensor imaging?<br />What metrics can be derived from DTI? |
| | 30 | [3.2 Constrained Spherical Deconvolution (CSD)][episode05] | What is Constrained Spherical Deconvolution (CSD)?<br />What does CSD offer compared to DTI? |
| **Tractography** | 30 | [4 Tractography][episode06] | What information can dMRI provide at the long range level? |
| | 30 | [4.1 Local tractography][episode07] | FIXME |
| | 30 | [4.1.1 Deterministic tractography][episode08] | FIXME |
| | 30 | [4.1.2 Probabilistic tractography][episode09] | Why do we need tractography algorithms beyond the deterministic ones?<br />How is probabilistic tractography different from deterministic tractography? |

## Contributing

We welcome all contributions to improve the lesson! Maintainers will do their best to help you if you have any
questions, concerns, or experience any difficulties along the way.

We'd like to ask you to familiarize yourself with our [Contribution Guide](CONTRIBUTING.md) and have a look at
the [more detailed guidelines][lesson-example] on proper formatting, ways to render the lesson locally, and even
how to write new episodes.

Please see the current list of [issues][link_issues] for ideas for contributing to this
repository. For making your contribution, we use the GitHub flow, which is
nicely explained in the chapter [Contributing to a Project](http://git-scm.com/book/en/v2/GitHub-Contributing-to-a-Project) in Pro Git
by Scott Chacon.
Look for the tag ![good_first_issue](https://img.shields.io/badge/-good%20first%20issue-gold.svg). This indicates that the maintainers will welcome a pull request fixing this issue.

## Maintainer(s)

Current maintainers of this lesson are

* [Jason Kai][jason_kai]
* [Olivia Stanley][olivia_stanley]
* [Michael Joseph][michael_joseph]
* [Jon Haitz Legarreta Gorroño][jon_legarreta]

## Authors

A list of contributors to the lesson can be found in [AUTHORS](AUTHORS)

## License

Instructional material from this lesson is made available under the Creative
Commons Attribution (CC BY 4.0) license. Except where otherwise noted, example
programs and software included as part of this lesson are made available under
the MIT license. For more information, see [LICENSE](LICENSE.md).

## Citation

To cite this lesson, please consult with [CITATION](CITATION)

[create_slack_svg]: https://img.shields.io/badge/Create_Slack_Account-The_Carpentries-071159.svg
[slack_heroku_invite]: https://swc-slack-invite.herokuapp.com
[slack_channel_status]: https://img.shields.io/badge/Slack_Channel-neuroimaging-E01563.svg
[slack_channel_url]: https://swcarpentry.slack.com/messages/CCJBHKCHZ
[binder_svg]: https://mybinder.org/badge_logo.svg
[binder_url]: https://mybinder.org/v2/gh/josephmje/SDC-BIDS-dMRI/gh-pages
[episode01]: https://carpentries-incubator.github.io/SDC-BIDS-dMRI/01-introduction_diffusion_data/index.html
[episode02]: https://carpentries-incubator.github.io/SDC-BIDS-dMRI/02-diffusion_preprocessing/index.html
[episode03]: https://carpentries-incubator.github.io/SDC-BIDS-dMRI/03-local_orientation_reconstruction/index.html
[episode04]: https://carpentries-incubator.github.io/SDC-BIDS-dMRI/04-diffusion_tensor_imaging/index.html
[episode05]: https://carpentries-incubator.github.io/SDC-BIDS-dMRI/05-constrained_spherical_deconvolution/index.html
[episode06]: https://carpentries-incubator.github.io/SDC-BIDS-dMRI/06-tractography/index.html
[episode07]: https://carpentries-incubator.github.io/SDC-BIDS-dMRI/07-local_tractography/index.html
[episode08]: https://carpentries-incubator.github.io/SDC-BIDS-dMRI/08-deterministic_tractography/index.html
[episode09]: https://carpentries-incubator.github.io/SDC-BIDS-dMRI/09-probabilistic_tractography/index.html
[lesson-example]: https://carpentries.github.io/lesson-example
[link_issues]: https://github.com/conp-pcno-training/SDC-BIDS-dMRI/issues
[jason_kai]: https://github.com/kaitj
[olivia_stanley]: https://github.com/ostanley
[michael_joseph]: https://github.com/josephmje
[jon_legarreta]: https://github.com/jhlegarreta
