# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt

from fury import window


def generate_axial_superior_slice_view(scene, size=(600, 600), offscreen=True):
    """Generate an axial superior view of the scene.

    Arguments
    ---------
    scene : window.Scene
        Scene to be transformed.
    size : Tuple (int, int)
        Size of each view.
    offscreen : bool
        `True` to prevent the scene window from being displayed.

    Returns
    -------
        A `window.Scene` instance seen from an axial superior view.
    """

    axial_scene = window.snapshot(scene, size=size, offscreen=offscreen)

    return axial_scene


def generate_sagittal_right_view(scene, size=(600, 600), offscreen=True):
    """Generate a sagittal right view of the scene.

    Arguments
    ---------
    scene : window.Scene
        Scene to be transformed.
    size : Tuple (int, int)
        Size of each view.
    offscreen : bool
        `True` to prevent the scene window from being displayed.

    Returns
    -------
        A `window.Scene` instance seen from a sagittal right view.
    """

    scene.yaw(90)
    scene.roll(-90)
    scene.reset_camera()
    sagittal_scene = window.snapshot(scene, size=size, offscreen=offscreen)

    # Reset scene
    scene.roll(90)
    scene.yaw(-90)

    return sagittal_scene


def generate_coronal_anterior_view(scene, size=(600, 600), offscreen=True):
    """Generate a coronal anterior view of the scene.

    Arguments
    ---------
    scene : window.Scene
        Scene to be transformed.
    size : Tuple (int, int)
        Size of each view.
    offscreen : bool
        `True` to prevent the scene window from being displayed.

    Returns
    -------
        A `window.Scene` instance seen from a coronal anterior view.
    """

    scene.yaw(90)
    scene.pitch(-90)
    scene.reset_camera()
    coronal_scene = window.snapshot(scene, size=size, offscreen=offscreen)

    # Reset scene
    scene.pitch(90)
    scene.yaw(-90)

    return coronal_scene


def generate_anatomical_volume_views(*actors, size=(600, 600)):
    """Generate anatomical (coronal anterior, sagittal right and axial
    superior) views of the actor(s).

    Arguments
    ---------
    actors : vtkActor
        Actor(s) to be displayed.
    size : Tuple (int, int)
        Size of each view.

    Returns
    -------
        A tuple of `window.Scene` instances containing the anatomical views of
        the actor(s).
    """

    offscreen = True

    # Create the scene
    scene = window.Scene()

    # Add the each data volume to the scene
    for actor in actors:
        scene.add(actor)

    axial_scene = generate_axial_superior_slice_view(
        scene, size=size, offscreen=offscreen)

    sagittal_scene = generate_sagittal_right_view(
        scene, size=size, offscreen=offscreen)

    coronal_scene = generate_coronal_anterior_view(
        scene, size=size, offscreen=offscreen)

    return coronal_scene, sagittal_scene, axial_scene


def generate_anatomical_slice_views(
        slices, *actors, size=(600, 600)):
    """Generate anatomical (coronal anterior, sagittal right and axial
    superior) views of the actor(s) at the provided slices.

    Arguments
    ---------
    slices: Tuple (int, int, int)
        Image data to be displayed.
    actors : vtkActor
        Actor(s) to be displayed.
    size : Tuple (int, int)
        Size of each view.

    Returns
    -------
        A tuple of `window.Scene` instances containing the anatomical views of
        the actor(s).
    """

    offscreen = True

    # Create the scene
    scene = window.Scene()

    # Add the each data volume to the scene
    for actor in actors:
        scene.add(actor)

    # Set the actor to the axial slice to be shown
    for actor in actors:
        actor.display(z=slices[2])

    axial_scene = generate_axial_superior_slice_view(
        scene, size=size, offscreen=offscreen)

    # Set the actor to the sagittal slice to be shown
    for actor in actors:
        actor.display(x=slices[0])

    sagittal_scene = generate_sagittal_right_view(
        scene, size=size, offscreen=offscreen)

    # Set the actor to the coronal slice to be shown
    for actor in actors:
        actor.display(y=slices[1])

    coronal_scene = generate_coronal_anterior_view(
        scene, size=size, offscreen=offscreen)

    return coronal_scene, sagittal_scene, axial_scene


def generate_figure(
        scene_titles, *scenes, cmap=None, figsize=(20, 20)):
    """Generate the figure containing the scenes.

    Arguments
    ---------
    scene_titles : list
        Scene title.
    *scenes : ndarray
        Scene data to be displayed.
    cmap : str
        Colormap to be applied to the overlay.
    figsize : Tuple (int, int)
        Subfigure size.

    Returns
    -------
        A `Figure` instance containing the provided scene(s).
    """

    origin = "lower"

    scene_count = len(scenes)

    fig, axes = plt.subplots(1, scene_count, figsize=figsize)

    # Plot the scenes
    for i, (scene, title) in enumerate(zip(scenes, scene_titles)):
        axes[i].imshow(scene, cmap=cmap, origin=origin)
        axes[i].axis("off")
        axes[i].set_title(title)

    return fig


def generate_anatomical_volume_figure(
        *actors, cmap=None, viewsize=(600, 600), figsize=(20, 20)):
    """Generate anatomical (coronal anterior, sagittal right and axial
    superior) views of the actor(s).

    Arguments
    ---------
    actors : vtkActor
        Actor(s) to be displayed.
    cmap : str
        Colormap to be applied to the overlay.
    viewsize : Tuple (int, int)
        Composite figure size.
    figsize : Tuple (int, int)
        Subfigure size.

    Returns
    -------
        A `Figure` instance containing the anatomical views of the actor(s).
    """

    scene_titles = ["Coronal anterior", "Sagittal right", "Axial superior"]

    coronal_scene, sagittal_scene, axial_scene = \
        generate_anatomical_volume_views(*actors, size=viewsize)

    fig = generate_figure(
        scene_titles, coronal_scene, sagittal_scene, axial_scene, cmap=cmap,
        figsize=figsize)

    return fig


def generate_anatomical_slice_figure(
        slices, *actors, cmap=None, viewsize=(600, 600), figsize=(20, 20)):
    """Generate anatomical (coronal anterior, sagittal right and axial
    superior) views of the actor(s) at the provided slices.

    Arguments
    ---------
    slices: Tuple (int, int, int)
        Image data to be displayed.
    actors : vtkActor
        Actor(s) to be displayed.
    cmap : str
        Colormap to be applied to the overlay.
    viewsize : Tuple (int, int)
        Composite figure size.
    figsize : Tuple (int, int)
        Subfigure size.

    Returns
    -------
        A `Figure` instance containing the anatomical views of the actor(s).
    """

    scene_titles = ["Coronal anterior", "Sagittal right", "Axial superior"]

    coronal_scene, sagittal_scene, axial_scene = \
        generate_anatomical_slice_views(slices, *actors, size=viewsize)

    fig = generate_figure(
        scene_titles, coronal_scene, sagittal_scene, axial_scene, cmap=cmap,
        figsize=figsize)

    return fig
