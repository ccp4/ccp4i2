#############
Xia2 with XDS
#############
   Xia2 is a software system that will perform automated X-ray diffraction
   data processing on your behalf. This particular pipeline will use XDS
   to index and integrate the reflection data.

   One important thing to note is that XDS is not part of the CCP4 software
   suite and will need to be installed indepedently. Details can be found on
   the XDS home-page. If you are a commercial user you will need an XDS
   license.

Setup
=====

   For the interface to work you will first needto install XDS, which is
   available for both Linux and Apple Mac from the XDS homepage (see below).
   Please follow the instructions there are test that XDS is working from
   the command line.

   After XDS is installed and working from the command line you will need
   to add the program to i2. You can do this by selecting 'edit - preferences'
   from the i2 tabs, and the selecting other programs. You should then add
   xds to the list, as shown below (note if you installed xds in a different
   folder you will need to adjust the path).

   |image1|

Input
=====

   In order to run the Xia2-Dials Interface you need to provide, at a
   minimum, a set of x-ray data images. DIALS supports a wide range of
   image formats and further details can be on the DIALS web-page.

   The image data can be selected for the task by either of two ways. The
   first is to select one individual image file (from a folder containing
   a set of images) & then specifying a range of images (eg. 1 to 100, if
   left blank all files are selected). The other way is to select the folder
   containing the images.

   Options for both Xia2 and XDS are also exposed

   -  `Xia2 Home Page <https://xia2.github.io/index.html>`__.
   -  `Xia2 Run Options <https://xia2.github.io/parameters.html>`__.
   -  `XDS Home Page <https://xds.mr.mpg.de/>`__.


.. |image1| image:: xia2xds_task_1.png
