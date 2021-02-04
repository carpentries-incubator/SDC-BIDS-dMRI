---
title: Setup
---

This lesson uses [Jupyter Notebook] as a web-based interactive computational
environment for learning. The `Jupyter Notebooks` can be run through [Binder],
which builds a computing environment with all of the necessary software
pre-installed. However, users may choose to run the notebooks locally. In that
case, there are several pieces of software that must be installed. Although the
main components required are `Python`-based, there are a few additional
non-`Python` tools required.

Instructors should allocate some time to have all components installed.

> ## Internet
> [Binder] requires users to have internet access to launch the `Jupyter`
> `Notebooks`. If users choose to run the notebooks locally, they will require
> internet access to download and install the software dependencies.
{: .prereq}

## Binder

[Binder] enables the user to open the collection of notebooks in this lesson in
a web-based executable and interactive environment. No additional software needs
to be installed locally in this case; it suffices to click on the
<kbd>launch binder</kbd> badge in the repository.

## Local

If users choose to run the `Jupyter Notebooks` locally, the following
dependencies will need to be installed:

- [ANTs] : used to register different anatomical data.
- [FSL] : used for different data preprocessing steps.
- [DIPY] : used for diffusion MRI data processing.
- [FURY] : used for anatomical data visualisation purposes.
- [Matplotlib] : used for data visualisation purposes.
- [Nilearn] : used for anatomical data visualisation purposes.
- [osfclient] : used to download the necessary data.
- [PyBIDS] : used to check the data structure [BIDS] compliance.

> ## Bash
>
> The installation instructions require terminal application (`bash`, `zsh`, or
> others) skills.
{: .prereq}

> ## Git
> In order to run the notebooks locally, users will need to retrieve the code
> using [Git] or directly as a zipped file from the lesson's [workshop-repo].
{: .prereq}

> ## Virtual environments
>
> A [virtual environment] is a tool that helps to keep dependencies required by
> different projects separate by creating an isolated, self-contained directory
> tree containing a set of particular `Python` packages required for a
> particular version of `Python`, and a project. This allows one to have
> different `Python` versions or different versions of a given package within
> the same computer without running the risk of mixing incompatible versions of
> different packages, and hence eliminating the risk of affecting the
> development or execution of a project.
>
> There are multiple package management tools (regardless of whether a virtual
> environment is used or not). These setup instructions will use [pip].
{: .callout}

For the setup purposes, it will be assumed that [Python] is already installed in
the local environment. Note that only Python 3+ versions are supported.

### Linux

`ANTs` and `FSL` are command line tools. Given the disparity and extent of the
steps involved in their installation, users should follow the specific
installation instructions in the corresponding official documentation pages.

The rest of the dependencies (`DIPY`, `FURY`, `Matplotlib`, `Nilearn`,
`osfclient`, and `PyBIDS`) of the lesson are `Python` packages. These also rely
on a number of other `Python` packages, which are automatically installed when
installing the former.

In order to install the `Python` dependencies, it suffices to run:
~~~
$ pip install -r requirements.txt
~~~
{: .bash}

from the root of the repository folder.

Users may choose to run the notebooks using [Jupyter] or [iPython]. The
`Jupyter` dependency can be installed by running:
~~~
$ pip install jupyter
~~~
{: .bash}

`iPython` can be installed through `pip` running:
~~~
$ pip install ipython
~~~
{: .bash}

Some `Python` distributions, such as the one provided by [Anaconda], might
ship these packages by default.

> ## Test the installation
>
> Test installation information for a package can be checked by running, for
> example:
> ~~~
> $ pip show dipy
> ~~~
> {: .bash}
>
> Similarly, it can be checked that a given package can be imported in Python by
> running, for example:
> ~~~
> $ python
> >>> import dipy
> ~~~
> {: .bash}
>
> Alternatively, the package version can also be checked by running, for example:
> ~~~
> $ python
> >>> import dipy
> >>> print(dipy.__version__)
> ~~~
> {: .bash}
>
> You can also see the packages and versions of all `pip`-installed dependencies
> by typing:
> ~~~
> $ pip freeze
> ~~~
> {: .bash}
{: .discussion}

In order to run the notebooks, the notebook server needs to be started. Once the
current directory changed to the root of the code directory, the server is
started running:
~~~
$ ipython notebook
~~~
{: .bash}

if using `IPython`, and running:

~~~
jupyter notebook
~~~
{: .bash}

if using `Jupyter`.

In either case, the commands will print some information about the notebook
server in the terminal, and a web browser will be opened to the URL of the web
application (by default, http://127.0.0.1:8888). The users will be presented to
the directory structure of the current directory, and they will be able to run
the notebook of interest.

For additional information about `Python` setups besides the package manuals,
users are encouraged to read the [Programming with Python] Carpentries lesson.

The data used in the lesson is hosted in [OSF]. It can be downloaded by running:
~~~
$ osf -p cmq8a clone ./data
~~~
{: .bash}

Notebooks expect them to be placed in the `data` folder that exists in the root
of the repository.


{% include links.md %}
