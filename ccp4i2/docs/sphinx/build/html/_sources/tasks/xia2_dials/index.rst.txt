###############
Xia2 with DIALS
###############
   Xia2 is a software system that will perform automated X-ray diffraction
   data processing on your behalf. This particular pipeline will use DIALS
   to index and integrate the reflection data.

Input
=====
   In order to run the Xia2-Dials Interface you need to provide, at a
   minimum, a set of x-ray data images. DIALS supports a wide range of
   image formats and further details can be on the DIALS web-page.

   |image1|

   The image data can be selected for the task by either of two ways. The
   first is to select one individual image file (from a folder containing
   a set of images) & then specifying a range of images (eg. 1 to 100, if
   left blank all files are selected). The other way is to select the folder
   containing the images.

   A large set of options for both Xia2 and DIALS are also made available for
   the user in the standard and advanced options tabs. These, for example,
   allow the user to manual set masks, and alter the algorithms DIALS uses
   for integration and spot finding.

   -  `Xia2 Home Page <https://xia2.github.io/index.html>`__.
   -  `Xia2 Run Options <https://xia2.github.io/parameters.html>`__.
   -  `DIALS Home Page <https://dials.diamond.ac.uk/>`__.

.. |image1| image:: xia2dials_task_1.png
