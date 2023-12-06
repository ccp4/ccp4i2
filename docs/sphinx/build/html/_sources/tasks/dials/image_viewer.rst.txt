==================
DIALS Image Viewer
==================

DIALS Image Viewer can be used for viewing diffraction images, optionally overlaid with the results of spot finding, indexing or integration. It is always worth inspecting the images with image viewer, especially if you are having issues with spot finding.

`DIALS JSON <https://dials.github.io/documentation/data_files.html#experiment-list-files>`_ and `Xia2 or DIALS pickle data file <https://dials.github.io/documentation/data_files.html#reflection-files>`_ need to be provided for this task.

Once the task is run, you will be directed to the main image viewer screen:

.. image:: image_1.png
         :scale: 50 %
         :align: center

Here you can adjust the brightness of images using the brightness slider (**1**) a little until you can clearly see the spots on the first image (something in the range of 10-20 should make the spots obvious). You can also change the color scheme, toggle various information markers like beam center, and try different configurations for the spot finding (**2**). You can use the image mode buttons (**3**) to preview how the parameters affect the spot-finding algorithm. The effect of changing the spot-finding parameters can be observed by displaying the threshold view of the image. This may be useful when commissioning a new type of detector or experiment.


**References**

`Winter, G., Waterman, D. G., Parkhurst, J. M., Brewster, A. S., Gildea, R. J., Gerstel, M., Fuentes-Montero, L., Vollmar, M., Michels-Clark, T., Young, I. D., Sauter, N. K. & Evans, G. (2018). Acta Cryst. D74, 85-97. <https://doi.org/10.1107/S2059798317017235>`_