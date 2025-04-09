##################
Pimple File Format
##################

Introduction
------------

The data file format of Pimple or QLoggraph is XML. The format is a
strict XML format which must pass `this
schema <http://www.ysbl.york.ac.uk/~mcnicholas/schema/CCP4ApplicationOutput.xsd>`__.
Everything is contained inside one element with the name
CCP4ApplicationOutput. Within the CCP4ApplicationOutput elements are
optionally any number of elements describing 2D plot data (CCP4Table),
3D contour plot data (Surface), and font preferences (Fonts). It is
possible that other preferences may be added in the future. The main
bulk of the file format definition is for the CCP4Table, which is an
XML-ified version of the $TABLE tables used in loggraph. The intention
is that, although the file format for Pimple is strict XML, it mirrors
very closely $TABLE such that changing programs to write Pimple XML
format should be a trivial change to CCP4 programs. The structure of the
tables and the headers is essentially the same, but simply wrapped up in
XML elements. The important changes are mostly in the description of
plots using such data.

Tags
----

-  CCP4ApplicationOutput (element)

   This is the outer XML element which contains all others. All pimple
   files should be of the form:

   ::

      <xml version="1.0" encoding="UTF-8" ?>
      <CCP4ApplicationOutput 
      xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
      xsi:noNamespaceSchemaLocation="http://www.ysbl.york.ac.uk/~mcnicholas/schema/CCP4ApplicationOutput.xsd">
      ⋮
      <CCP4ApplicationOutput>

   -  CCP4Table (element)

      This element is the analogue of loggraph's $TABLE. A Pimple input
      file may contain any number of this element, each instance will
      define a set of data to plot (as selectable by the left hand
      combo-box in Pimple). Each CCP4Table may contain any number of
      plot elements, each of which describes a single graph with any
      number of plotlines. The plot elements provide entries in the
      right-hand combo-box in Pimple.

      This element may have the following:

      -  title (attribute, string, optional). The name of the table as
         it will appear in the left-hand combo-box.
      -  headers (element, string, optional) The column headers for the
         data, i.e. a label for each column of the data. This is the
         equivalent of the first $$…$$ stanza in a $TABLE. The header
         labels may be separated by spaces or any other character as
         given in the following attribute.

         -  separator (attribute, string, optional)

      -  data (element, string, optional) The data, like the headers,
         the data may be separated by a space or a character specified
         by the separator attribute. The format is very similar to that
         of the second $$…$$ stanza of $TABLE, except that columns may
         not run into each other, data items must be separated.

         -  separator (attribute, string, optional)

      -  plot (element) A graph or plot, the intended final picture.
         This consists of any number of plotlines, lines, circles and
         polygons. This is the most complex element in the Pimple
         format, with many elements describing the appearance of the
         plot. All of the following elements may occur any number of
         times, including 0.

         -  plottype (element, string) The type of plot, currently only
            xy is sensible.

         -  title (element, string) The title as it will appear at top
            of plot and in right-hand combo-box.

         -  xlabel (element, string) The x-axis label

         -  ylabel (element, string) The y-axis label

         -  rylabel (element, string) The y-axis label on a right-hand
            axis.

         -  xscale (element, string) A scaling may be applied to the
            data in the x-dirtection. The only currently supported value
            is oneoversqrt. Pimple usually understands the column
            headers from CCP4 programs and so this element is not
            generally necessary.

         -  yscale (element, string) See xscale above.

         -  xrange (element) A range of values for the extents of the
            x-axis, 2 compulsory attributes:

            -  min (attribute, float, mandatory)
            -  max (attribute, float, mandatory)

         -  yrange (element) See xrange above.

            -  rightaxis (attribute, boolean) If true, the range is for
               a second right-hand axis on the graph. This give a
               "dual-axis" appearence.

         -  xbreaks (element) The axes may be broken in an arbitrary
            number of places. It is (currently) only sensible to break
            either x or y axis, not both. A break is defined by a range
            in a break element:

            -  break (element) The range of the break, 2 compulsory
               attributes:

               -  min (attribute, float)
               -  max (attribute, float)

         -  ybreaks (element) See xbreaks above.

         -  legendposition (element) The position of the legend in
            fractional screen coordinates: (0,0) is upper left; (1,1) is
            lower right. 2 compulsory attributes:

            -  x (attribute, float, mandatory)
            -  y (attribute, float, mandatory)

            If this element is not present, Pimple will try to place the
            legend somewhere sensible.

         -  showlegend (element, boolean) Whether or not to show the
            legend.

         -  xintegral (element, boolean) Whether or not the x-axis data
            are to be treated as integers. Default is false and all data
            are treated as floats.

         -  yintegral (element, boolean) See xintegral above.

         -  yintegral (element, boolean) See xintegral above.

         -  plotline (element)

            An individual plot line representing one column of the data
            plotted against another. Which columns of data correspond to
            x and which to y are given by the compulsory xcol and ycol
            attributes.

            -  xcol (attribute, integer, mandatory) Which column of the
               data element is to be used for x values.
            -  ycol (attribute, integer, mandatory) Which column of the
               data element is to be used for y values.
            -  symbol (element,string) The symbol which will be plotted
               at each data point for this line. Allowed values are:
               'o','s','d','^','>','<','p','h','*','$\lambda$',
               '$\bowtie$', '$\circlearrowleft$', '$\clubsuit$',
               '$\checkmark$'. These values correspond to:
               'Circle','Square','Diamond','Up arrow','Right
               arrow','Left arrow','Pentagon','Hexagon','Star',
               'Lambda','Bow tie','Circle arrow left','Clubs (as in
               playing cards)','Tick'. See matplotlib documentation for
               more information. If plotline elements do not specify
               symbol then the symbols just described will be used for
               each plotline in order.
            -  symbolsize (element, integer) The size of the symbol.
            -  linestyle (element, string) The style of the line,
               allowed values: '-','--','-.',':','.', corresponding to:
               'Solid','Dashed','Dash-dot','Dotted','Blank'. See
               matplotlib documentation for further details.
            -  linesize (element, integer) The width of the line.
            -  colour (element, string) The colour of the line and
               symbol. Allowed values red, green, blue, yellow, magenta,
               cyan, black, or r, g, b, y, m, c, k, or any other colour
               specification understood by matplotlib, e.g. orange,
               #ff7700.
            -  label (element, string) The label that will be used for
               this line in the legend.
            -  visible (element, boolean) Whether or not to show this
               line.
            -  showinlegend (element, boolean) Whether or not the line
               will appear in the legend.
            -  rightaxis (attribute, boolean) If true, the y-coordinate
               of the plot corresponds to the range on the second
               right-hand axis of the graph.

         -  circle (element) The definition of a circle to place onto
            the picture.

            -  radius (attribute, madatory) The radius of the circle in
               data units.
            -  xpos (attribute, madatory) The x position of the circle
               in data units.
            -  ypos attribute, madatory) The y position of the circle in
               data units.
            -  linestyle (attribute, string) See plotline->linestyle.
            -  linesize (attribute, integer) See plotline->linesize.
            -  linecolour (attribute, string) The colour of the outline.
               See plotline->colour.
            -  fillcolour (attribute, string) The colour which fills the
               circle, if absent the circle will not be filled. See
               plotline->colour.

         -  line (element) A simple line to draw on the picture, not
            associated with any data.

            -  x1 (attribute, madatory) The start x position of the
               line.
            -  x2 (attribute, madatory) The start y position of the
               line.
            -  y1 (attribute, madatory) The end x position of the line.
            -  y2 (attribute, madatory) The end y position of the line.
            -  linestyle (attribute, string) See plotline->linestyle.
            -  linesize (attribute, integer) See plotline->linesize.
            -  linecolour (attribute, string)See plotline->colour.

         -  polygon (element) A polygon to draw on the picture. The
            coordinates of the vertices are of the form x1, y1, x2, y2,
            x3, y3, x4, y4,…, xn, yn. These are data coordinates. The
            format of this data is almost certain to change in the
            future.

            -  linestyle (attribute, string) See plotline->linestyle.
            -  linesize (attribute, integer) See plotline->linesize.
            -  linecolour (attribute, string) The colour of the outline.
               See plotline->colour.
            -  fillcolour (attribute, string) The colour which fills the
               polygon, if absent the polygon will not be filled. See
               plotline->colour.

      A complete example. This example has just one plot (option in
      Pimple's right-hand combo-box). The plot has 4 plot lines, with
      custom colour, style, symbol, etc. There are 3 x-axis breaks; a
      range is specified for the y-axis; the xscaling is 1/√ x ; &#xC5;
      is the HTML entity for Å.

      ::

         <CCP4Table title="Cycle    1. Rfactor analysis, F distribution v resln">
         <headers separator=" "> 
         M(4SSQ/LL) NR_used %_obs M(Fo_used) M(Fc_used) Rf_used WR_used
         NR_free M(Fo_free) M(Fc_free) Rf_free   WR_free
         </headers>
         <data>
          0.019     723  84.15   292.1   292.4  0.21  0.22      52   293.1   299.2  0.19  0.19
          0.057    1454  98.70   350.9   346.1  0.17  0.18      66   332.3   329.1  0.14  0.14
          0.095    1807  97.86   276.4   266.8  0.20  0.19     113   293.0   284.9  0.16  0.15
          0.132    2143  97.91   188.7   185.7  0.24  0.27     103   177.6   178.5  0.26  0.28
          0.170    2335  97.14   155.2   150.3  0.25  0.26     144   167.7   154.0  0.21  0.23
          0.208    2577  96.12   140.9   129.5  0.26  0.27     151   138.2   130.8  0.27  0.28
          0.245    2809  96.76   119.6   107.5  0.28  0.28     151   126.5   110.4  0.26  0.27
          0.283    2916  93.73    93.9    87.9  0.30  0.29     165    87.9    84.5  0.34  0.31
          0.321    2931  88.75    76.8    72.5  0.32  0.31     155    75.8    70.6  0.31  0.28
          0.358    2874  83.09    60.8    59.3  0.34  0.32     162    62.5    58.9  0.33  0.32
          0.396    3022  82.69    54.0    52.2  0.34  0.32     173    53.6    53.3  0.37  0.34
          0.434    3170  83.43    48.0    45.8  0.35  0.31     178    44.7    44.0  0.39  0.33
          0.471    3248  81.57    42.4    40.8  0.36  0.32     173    41.2    39.0  0.34  0.30
          0.509    3316  80.92    38.2    36.2  0.37  0.32     192    38.4    36.9  0.34  0.29
          0.547    3417  79.43    34.6    33.0  0.38  0.33     150    35.5    33.6  0.37  0.34
          0.584    3427  77.12    32.3    30.1  0.40  0.35     166    33.7    30.7  0.41  0.36
          0.622    3494  76.72    30.4    28.0  0.39  0.35     170    29.9    27.2  0.36  0.31
          0.660    3493  74.92    29.1    25.7  0.39  0.36     187    28.7    26.0  0.37  0.33
          0.697    3544  74.44    27.3    23.9  0.40  0.37     205    26.1    24.3  0.37  0.35
          0.735    3401  69.25    26.3    22.5  0.41  0.38     182    23.2    20.7  0.43  0.41
         </data>
         <plot>
         <xscale>oneoversqrt</xscale>
         <plottype>xy</plottype>
         <title>
         Cycle    1. <Rfactor> v. resln
         </title>
         <xlabel>
         Resolution / &#xC5;
         </xlabel>
         <ylabel>
         R factor
         </ylabel>
         <yrange min="-.1" max="0.7"/>
         <legendposition x="0.75" y="0.05"/>
         <xbreaks>
         <break min="0.24" max="0.26"/>
         <break min="0.3" max="0.35"/>
         <break min="0.4" max="0.5"/>
         </xbreaks>
         <plotline xcol="1" ycol="6">
         <colour>red</colour>
         <linesize>2</linesize>
         </plotline>
         <plotline xcol="1" ycol="7">
         <colour>blue</colour>
         <linesize>1</linesize>
         </plotline>
         <plotline xcol="1" ycol="11">
         <colour>red</colour>
         <linesize>3</linesize>
         </plotline>
         <plotline xcol="1" ycol="12">
         <colour>blue</colour>
         <linesize>3</linesize>
         <symbol>*</symbol>
         <symbolsize>10</symbolsize>
         </plotline>
         </plot>

   -  Surface (element)

      This element defines a 3D contour plot. The various style, colour,
      etc. options which are available in 2D xy plots are not in this
      element: matplotlib is more or less left to to present the data in
      its default manner. The data are presented as x1, y1, z1, x2, y2,
      z2, x3, y3, z3, x4, y4, z4, …, xn, yn, zn.

      -  rows (attribute, integer, mandatory) The number of rows.
      -  column (attribute, integer, mandatory) The number of columns.
      -  title (attribute, string) The title.

      ::

   -  Fonts (element)

      This element represents options for the fonts.

      -  titleFont (element) The font used for the titles at the top of
         plots.

         -  family (attribute)
         -  size (attribute)
         -  weight (attribute)
         -  slant (attribute)

      -  legendFont (element) The legend font
      -  axesTickerFont (element) The font used for the numbers on the
         axes tickers.
      -  axesLabelFont (element) The font used for the axes labels.

      Example:

      ::

         <Fonts>
         <titleFont family="Bitstream Vera Sans" size="20"/>
         <Fonts>

A complete example file
-----------------------

::

   <?xml version="1.0" encoding="UTF-8" ?>
   <CCP4ApplicationOutput 
   xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
   xsi:noNamespaceSchemaLocation="http://www.ysbl.york.ac.uk/~mcnicholas/schema/CCP4ApplicationOutput.xsd">
   <Fonts>
   <titleFont family="Bitstream Vera Sans" size="20"/>
   </Fonts>

   <CCP4Table title="Rfactor analysis, stats vs cycle">
   <headers separator=" "> 
       Ncyc    Rfact    Rfree     FOM      -LL     -LLfree  rmsBOND  zBOND rmsANGL  zANGL rmsCHIRAL
   </headers>
   <plot>
   <title>wibble</title>
   <xintegral>true</xintegral>
   <xrange min="0" max="6"/>
   <plotline xcol="1" ycol="2">
   <label>hmm</label>
   <colour>#ff7700</colour>
   </plotline>
   <plotline xcol="1" ycol="3">
   </plotline>
   <xlabel>
   Ncyc
   </xlabel>
   </plot>
   <data>
          0   0.2783   0.2668   0.741    250301.   13484.0   0.0284  1.426   4.570  2.469   0.170
          1   0.2270   0.2359   0.805    241026.   13136.4   0.0299  1.371   2.898  1.470   0.163
          2   0.2057   0.2194   0.837    235671.   12917.0   0.0307  1.373   2.642  1.330   0.169
          3   0.1950   0.2111   0.854    232521.   12787.8   0.0323  1.430   2.562  1.296   0.173
          4   0.1898   0.2073   0.862    230729.   12713.8   0.0327  1.434   2.485  1.259   0.175
          5   0.1873   0.2056   0.867    229749.   12676.0   0.0325  1.413   2.409  1.221   0.175
          6   0.1861   0.2048   0.870    229236.   12657.3   0.0320  1.380   2.338  1.184   0.174
          7   0.1854   0.2044   0.871    228920.   12648.3   0.0321  1.380   2.312  1.172   0.177
          8   0.1849   0.2044   0.872    228742.   12645.4   0.0322  1.384   2.301  1.168   0.181
          9   0.1848   0.2045   0.873    228632.   12645.6   0.0323  1.388   2.296  1.166   0.187
         10   0.1849   0.2045   0.873    228629.   12648.1   0.0317  1.352   2.254  1.142   0.189
   </data>
   </CCP4Table>

   <CCP4Table title="Cycle    1. Rfactor analysis, F distribution v resln">
   <headers separator=" "> 
   M(4SSQ/LL) NR_used %_obs M(Fo_used) M(Fc_used) Rf_used WR_used
   NR_free M(Fo_free) M(Fc_free) Rf_free   WR_free
   </headers>
   <data>
    0.019     723  84.15   292.1   292.4  0.21  0.22      52   293.1   299.2  0.19  0.19
    0.057    1454  98.70   350.9   346.1  0.17  0.18      66   332.3   329.1  0.14  0.14
    0.095    1807  97.86   276.4   266.8  0.20  0.19     113   293.0   284.9  0.16  0.15
    0.132    2143  97.91   188.7   185.7  0.24  0.27     103   177.6   178.5  0.26  0.28
    0.170    2335  97.14   155.2   150.3  0.25  0.26     144   167.7   154.0  0.21  0.23
    0.208    2577  96.12   140.9   129.5  0.26  0.27     151   138.2   130.8  0.27  0.28
    0.245    2809  96.76   119.6   107.5  0.28  0.28     151   126.5   110.4  0.26  0.27
    0.283    2916  93.73    93.9    87.9  0.30  0.29     165    87.9    84.5  0.34  0.31
    0.321    2931  88.75    76.8    72.5  0.32  0.31     155    75.8    70.6  0.31  0.28
    0.358    2874  83.09    60.8    59.3  0.34  0.32     162    62.5    58.9  0.33  0.32
    0.396    3022  82.69    54.0    52.2  0.34  0.32     173    53.6    53.3  0.37  0.34
    0.434    3170  83.43    48.0    45.8  0.35  0.31     178    44.7    44.0  0.39  0.33
    0.471    3248  81.57    42.4    40.8  0.36  0.32     173    41.2    39.0  0.34  0.30
    0.509    3316  80.92    38.2    36.2  0.37  0.32     192    38.4    36.9  0.34  0.29
    0.547    3417  79.43    34.6    33.0  0.38  0.33     150    35.5    33.6  0.37  0.34
    0.584    3427  77.12    32.3    30.1  0.40  0.35     166    33.7    30.7  0.41  0.36
    0.622    3494  76.72    30.4    28.0  0.39  0.35     170    29.9    27.2  0.36  0.31
    0.660    3493  74.92    29.1    25.7  0.39  0.36     187    28.7    26.0  0.37  0.33
    0.697    3544  74.44    27.3    23.9  0.40  0.37     205    26.1    24.3  0.37  0.35
    0.735    3401  69.25    26.3    22.5  0.41  0.38     182    23.2    20.7  0.43  0.41
   </data>
   <plot>
   <xscale>oneoversqrt</xscale>
   <plottype>xy</plottype>
   <title>
   Cycle    1. <Rfactor> v. resln
   </title>
   <xlabel>
   Resolution / &#xC5;
   </xlabel>
   <ylabel>
   R factor
   </ylabel>
   <yrange min="-.1" max="0.7"/>
   <legendposition x="0.75" y="0.05"/>
   <xbreaks>
   <break min="0.24" max="0.26"/>
   <break min="0.3" max="0.35"/>
   <break min="0.4" max="0.5"/>
   </xbreaks>
   <plotline xcol="1" ycol="6">
   <colour>red</colour>
   <linesize>2</linesize>
   </plotline>
   <plotline xcol="1" ycol="7">
   <colour>blue</colour>
   <linesize>1</linesize>
   </plotline>
   <plotline xcol="1" ycol="11">
   <colour>red</colour>
   <linesize>3</linesize>
   </plotline>
   <plotline xcol="1" ycol="12">
   <colour>blue</colour>
   <linesize>3</linesize>
   <symbol>*</symbol>
   <symbolsize>10</symbolsize>
   </plotline>
   </plot>

   <plot>
   <xscale>oneoversqrt</xscale>
   <showlegend>false</showlegend>
   <plottype>xy</plottype>
   <title>
   Cycle    1. <Fobs> and <Fc> v. resln
   </title>
   <xlabel>
   Resolution / &#xC5;
   </xlabel>
   <legendposition x="0.75" y="0.75"/>
   <plotline xcol="1" ycol="4">
   <colour>orange</colour>
   <linesize>2</linesize>
   <symbol>*</symbol>
   <symbolsize>20</symbolsize>
   </plotline>
   <plotline xcol="1" ycol="5">
   <colour>orange</colour>
   <linesize>2</linesize>
   <symbol>*</symbol>
   <symbolsize>20</symbolsize>
   </plotline>
   <plotline xcol="1" ycol="9">
   <colour>orange</colour>
   <linesize>2</linesize>
   <symbol>*</symbol>
   <symbolsize>20</symbolsize>
   </plotline>
   <plotline xcol="1" ycol="10">
   <colour>orange</colour>
   <linesize>2</linesize>
   <symbol>*</symbol>
   <symbolsize>20</symbolsize>
   </plotline>
   </plot>

   <plot>
   <xscale>oneoversqrt</xscale>
   <showlegend>false</showlegend>
   <plottype>xy</plottype>
   <title>
   Cycle    1. % observed v. resln
   </title>
   <xlabel>
   Resolution / &#xC5;
   </xlabel>
   <legendposition x="0.75" y="0.75"/>
   <plotline xcol="1" ycol="3">
   <colour>red</colour>
   <linesize>2</linesize>
   <linestyle>-.</linestyle>
   </plotline>
   </plot>


   </CCP4Table>
   <CCP4Table title="Cycle    1. Fom(<cos(DelPhi)>acentric, centric, overall v resln">
   <headers> 
    <4SSQ/LL> NREFa  FOMa  NREFc FOMc NREFall FOMall  SigmaA_Fc1
   </headers>
   <data separator=" ">
     0.0192   518   0.919   205   0.855   723   0.901  0.981
     0.0569  1173   0.949   281   0.872  1454   0.935  0.980
     0.0946  1535   0.922   272   0.878  1807   0.915  0.980
     0.1323  1872   0.863   271   0.812  2143   0.857  0.974
     0.1699  2078   0.847   257   0.779  2335   0.839  0.964
     0.2076  2318   0.852   259   0.825  2577   0.850  0.991
     0.2453  2543   0.814   266   0.748  2809   0.807  0.996
     0.2830  2677   0.790   239   0.748  2916   0.787  0.963
     0.3206  2720   0.760   211   0.690  2931   0.755  0.911
     0.3583  2692   0.730   182   0.630  2874   0.724  0.868
     0.3960  2843   0.705   179   0.663  3022   0.703  0.845
     0.4337  2988   0.724   182   0.646  3170   0.719  0.853
     0.4713  3081   0.743   167   0.698  3248   0.741  0.853
     0.5090  3146   0.701   170   0.661  3316   0.699  0.836
     0.5467  3256   0.691   161   0.670  3417   0.690  0.826
     0.5844  3290   0.623   137   0.600  3427   0.623  0.800
     0.6220  3345   0.691   149   0.661  3494   0.689  0.835
     0.6597  3355   0.678   138   0.624  3493   0.676  0.839
     0.6974  3402   0.688   142   0.641  3544   0.686  0.853
     0.7351  3270   0.662   131   0.616  3401   0.661  0.836
   </data>

   <plot>
   <yrange min="-.6" max="1"/>
   <xrange min="-1" max="10"/>
   <plottype>xy</plottype>
   <title>
   Cycle    1. % observed v. resln
   </title>
   <plotline xcol="1" ycol="3"/>
   <plotline xcol="1" ycol="5"/>
   <plotline xcol="1" ycol="7"/>
   <plotline xcol="1" ycol="8"/>
   <line x1="0" y1="0" x2="1" y2="1" linesize="5" linestyle="dashed"/>
   </plot>

   </CCP4Table>
   </CCP4ApplicationOutput>
