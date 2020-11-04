# Introduction to MRI and BIDS

![Build Status](https://github.com/carpentries-incubator/SDC-BIDS-dMRI/actions?query=workflow%3A"Build%2C+test")
[![Create a Slack Account with us][create_slack_svg]][slack_heroku_invite]
[![Slack Status][slack_channel_status]][slack_channel_url]
[![Binder][binder_svg]][binder_url]

An introduction to magnetic resonance imaging analysis in Python.

## Why Python?

Python is rapidly becoming the standard language for data analysis, visualization and automated workflow building. It is a free and open-source software that is relatively easy to pick up by new programmers. In addition, with Python packages such as `Jupyter` one can keep an interactive code journal of analysis - this is what we'll be using in the workshop. Using Jupyter notebooks allows you to keep a record of all the steps in your analysis, enabling transparency and ease of code sharing.

Another advantage of Python is that it is maintained by a large user-base. Anyone can easily make their own Python packages for others to use. Therefore, there exists a *very* large codebase for you to take advantage of for your neuroimaging analysis; from basic statistical analysis, to brain visualization tools, to advanced machine learning and multivariate methods!

## About the Lesson

This lesson teaches:
- a (re?) introduction to MR nomenclature - with BIDS
- "converting" your data to BIDS
- BIDS apps
- queueing up neuroimaging pipelines
- how neuroimaging data is stored

## Episodes

| # |  Episode | Time | Question(s) |
|--:|:---------|:----:|:------------|
| 1 | [Introduction to Diffusion and the Dataset][episode01] | 30 | How can I do the same operations on many different values? |
| 2 | [Diffusion Preprocessing][episode02] | 30 | How is MRI data organized |
| 3 | [Diffusion Tensor Imaging][episode03] | 30 | How can I do the same operations on many different files? |
| 4 | [Tractography][episode04] | 30 | How can I store many values together? |

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
[episode01]: https://conp-pcno-training.github.io/SDC-BIDS-dMRI/01-neuroimaging-fundamentals/index.html
[episode02]: https://conp-pcno-training.github.io/SDC-BIDS-dMRI/02-anatomy-of-nifti/index.html
[episode03]: https://conp-pcno-training.github.io/SDC-BIDS-dMRI/03-brain-imaging-data-structure/index.html
[episode04]: https://conp-pcno-training.github.io/SDC-BIDS-dMRI/04-open-mri-datasets/index.html
[lesson-example]: https://carpentries.github.io/lesson-example
[link_issues]: https://github.com/conp-pcno-training/SDC-BIDS-dMRI/issues
[jason_kai]: https://github.com/kaitj
[olivia_stanley]: https://github.com/ostanley
[michael_joseph]: https://github.com/josephmje
