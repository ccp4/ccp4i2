############################
CCP4i2 — CCP4i2 Search Tools
############################

**Please note, the CCP4i2 search tools are being radically changed. This
documentation is likely no longer valid.**

Job search
~~~~~~~~~~

The *Job search* tool on the *Project viewer* toolbar will find a job in
the project based on its name, when it was run and text in its title or
annotation (see `text search syntax <#text_search_syntax>`__. In
addition you can search:

-  the *Data history* for which you must select one job and the data
   type of input or output data - this will select preceeding jobs that
   provided input data and succeeding that used output data. Note that
   the job context menu (right mouse click on a job) in the *Job list*
   has a *Data history* option with similar functionality.
-  the jobs that used specific *control parameters* in running a job
   (beware this searches the parameter files which can be slow). See
   `details <#control_params>`__.
-  where an external *Imported file* is first used - you must select or
   provide the path name of the file

The 'hit' jobs are highlighted on the *Job list* in red. The *Data
history* search may well find more than one path after the job and you
can scroll between the options with the slider that is displayed above
the *Job list*.

Text search syntax
~~~~~~~~~~~~~~~~~~

The job search tool will search the job title (as it appears in the *Job
list*) or *Comments* (as entered in the *Comments* tab for the job) for
the text string which may contain the following wildcards:

+---------------------------+-----------------------------------------+
| Wildcard                  | Description                             |
+===========================+=========================================+
| %                         | A substitute for zero or more           |
|                           | characters                              |
+---------------------------+-----------------------------------------+
| \_                        | A substitute for a single character     |
+---------------------------+-----------------------------------------+
| [charlist]                | Sets and ranges of characters to match  |
+---------------------------+-----------------------------------------+
| [^charlist] or[!charlist] | Matches only a character NOT specified  |
|                           | within the brackets                     |
+---------------------------+-----------------------------------------+

For example to find jobs whose title begin with 'Refinement' you will
need to enter *Refinement%*

Searching control parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This tool is intended to help you find jobs that were run with specific
parameters. The job parameters are not saved in the database and so the
search must open the *input_params.xml* file for each searched job -
this may be slow. In the search dialog the taskname must be chosen and
the *Control parameters* radio button clicked to show the options. The
parameter name must be selected from a list - if you do not know the
parameter name open a job of the required task and hover the cursor over
the input field you are interested in; the tooltip that appears after a
short time should have the name of the parameter at the end of the help
text. The parameter name is usually uppercase letters. The definition of
the required value can have one or more components which are separated
by spaces so for example to find a value between 10 and 20 the
definition would be '*>10 <20*'. The following comparison operators can
be used:

========= =====================
Operation Description
========= =====================
<         less than
>         greater than
=         equal to
!=        not equal to
<=        less than or equal
>=        greater than or equal
========= =====================

If there is no operator then it is assumed to be 'equal to'. If a
comma-separated list of values is entered then any of those values will
be found.

Search projects
~~~~~~~~~~~~~~~

There is a simple text search tool at the bottom of the *Manage
projects* window. By default the text search is applied to the project
name, description and tags (the last two can be entered when creating
the project or using the *Edit description* tool in the *Manage
projects* window) but unchecking the appropriate boxes can remove any of
these from the search. The `syntax for the search
text <#text_search_syntax>`__ is the same as for the job search tool.
Hit projects are highlighted in red in the project list and any 'parent'
project with hit children is highlighted in pink.

| Click on the *Advanced search tools* button |Advanced search button|
  in the *Search* bar to see more search options. The search can be
  restricted to:
| projects active in a given period - the 'active' period of a project
  is from when it was first created until when it was last opened
| sub-projects of a specified project.
| It is also possible to find the projects that have imported a
  specified file or use a specified project directory.

Last modified: Wed Nov 16 17:44:04 GMT 2016

.. |CCP4 web logo| image:: ccp4i2.png
   :width: 75px

.. |Advanced search button| image:: ../../../../qticons/search-plus.png
   :width: 30px
   :height: 30px
