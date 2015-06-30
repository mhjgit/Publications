% $Id: README.txt,v 1.16 2005/12/01 16:16:27 frank Exp $

This distribution implements a class for the conference proceedings
which have layouts close to those of the AIP proceedings.

Currently it supports:

  AIP Conference Proceedings  6in   x 9in single column 
  AIP Conference Proceedings  8.5in x 9in single column 
  AIP Conference Proceedings  8.5in x 9in double column 
  ARLO --- Acoustics Research Letters Online

The features of the class are described in aipguide (.pdf or .tex).
Special requirements or restrictions for the above layouts are given
at the end of this guide!


Templates:
==========

The distribution contains one or more template files to be used as a
model for writing new articles. The are named template-XX.tex where
(XX depends on the proceeding layout, e.g. 6s for 6x9inch single).



General requirements and restrictions:
======================================

This class was designed to work with LaTeX release 1999/6/01 or a
later version. However releases from late 1996 onwards will work
correctly if the installation is complete, ie contains all required
files which may not be the case.

Releases earlier than 1996/12/01 might compile the documents but will
*not* result in correct output due to some bugs in these earlier
releases! (In addition you need to have at least a newer version of
docstrip.tex to install the aipproc distribution.)

With the exception of the packages natbib.sty and url.sty the class
only requires files which are part of a standard LaTeX distribution,
i.e., it will work if your installation contains the following
standard components: BASE, TOOLS, GRAPHICS, PSNFSS and the above
packages.

If your installation does not fulfill the above requirements you may
have to upgrade your LaTeX to be able to use this class.

You can test your installation by running the file aipcheck.tex
through your version of LaTeX. It will try to determine if everything
necessary is available and if not will make recommendations what can
be done about it. In certain cases you might be able to use the class
if you follow the suggestions, in other cases the only solution is to
upgrade your LaTeX.

The most recent LaTeX as well as natbib.sty and url.sty can be
obtained from CTAN sites (Comprehensive TeX Archive Network).

Refer to http://www.tug.org for more information on CTAN and
TeX in general.

A ready to run TeX system for various platforms which has
everything required is available on CD-ROM, look into
http://www.tug.org/texlive.html.

For loading individual packages from a CTAN archive refer to
http://www.ctan.org and search for the package name. Please omit
extensions such as .sty when searching, e.g., search for natbib rather
than natbib.sty, as such packages are often distributed in source form
only, e.g., as a .dtx file.

It is also possible to download a complete TeX/LaTeX installation from
CTAN, e.g., Miktex + Winedit + Ghostview. Finally, it is also possible
to download a CD-ROM image of the TeX-live CD from CTAN (roughly
300MB): search for 'texlive' (and make sure you select a suitable
mirror near you).


Class Documentation:
====================

The class documentation is provided as a PDF file but you can also
produce it from its source.

To produce the general guide, run aipguide.tex three times through LaTeX.
                                               ^^^^^^^^^^^
To produce the code documentation run the aipproc.dtx or the
fix2col.dtx file through LaTeX.



Files in the distribution:
==========================

Class file:
-----------

aipproc.cls	The main class file

aipxfm.sty	Front matter part of aipproc

aip-*.clo	Several layout files that implement layout styles for
		different proceeding styles


Included package files used by class:
-------------------------------------

fix2col.sty	This package is only included for people with older
		LaTeX installations and is not necessary with an
		up-to-date release.


BibTeX support:
---------------

aipproc.bst	BibTeX style file implementing something close the 
		requested layout for the AIP proceedings --- 
                to be used if natbib is available.


aipprocl.bst	BibTeX style file implementing something close the 
		requested layout for the AIP proceedings ---
		should only be used if natbib is not installed!




Templates:
----------

template-xx.tex Special template file which can be used as a model for
		writing a new article. The "xx" depends on the
		proceeding layout and might be 6s or arlo... 



Documentation files:
--------------------

Misc files:
-----------

README.txt	this file

aipguide.tex	Guide explaining all the features of the aipproc class.
		To successfully process this file a reasonable new LaTeX
		release is required (see section "General requirements and 
		restrictions" below).

aipguide-generic.pdf	Processed version of aipguide in PDF format. You
		should be able to view this with Adobe's Acrobat
		Reader which is available free of charge for most
		platforms.

aipguide-xx.pdf Special versions of the aipguide related to a specific
		proceeding layout, e.g., aipguide-arlo.pdf would describe
		only features needed the ARLO proceeding.
		These guides are not part of the distribution but
		might be provided.

FAQ.txt		Frequently Asked Questions

ChangeLog.fmi	Documentation of changes to the files in the distribution.





Documented sources and install scripts:
---------------------------------------

aipproc.dtx	Documented source of the main class file
aipparms.dtx	Documented source of the parameter settings

fix2col.dtx	Source version v0.03 of this package.

aipdoc.cls	Used for code documentation.

aipproc.ins	file to unpack aipproc.dtx and fix2col.dtx
aipparms.ins	file to unpack the .clo files holding the parameter settings



Test files:
-----------

aipcheck.tex	Checks the LaTeX installation for possible problems,
		used by sample.tex and aipguide.tex.

sample.tex 	A sample article that can be used to produce the
		various layouts supported by this class.

sample.bib	bib file for sample.tex

*.eps		a few PS files used as figures in sample.tex and aipguide.tex





Major packages used by class but not included:
----------------------------------------------

calc.sty	part of standard LaTeX (tools distribution)

graphicx.sty	part of standard LaTeX (graphics distribution)

mathptm.sty	implements Times as math fonts. Part of standard LaTeX
mathptmx.sty    (psnfss distribution)

mathtime.sty	implements commercial MathTime fonts as math fonts.
		package is available on CTAN, fonts must be purchased
		from Y&Y. (Package only used when requested by option.)

fontenc.sty	part of standard LaTeX

times.sty	part of standard LaTeX (psnfss distribution)

textcomp.sty	part of standard LaTeX

natbib.sty	implements extended citation. Should be available on
		any good LaTeX installation.

url.sty		implements URL support. Should be available on
		any good LaTeX installation.


Installation:
=============

Check if you have a later version of fix2col.sty, if so delete the
version coming with this distribution.

You may have to move or copy all files with the extensions .cls .clo
and .sty into a directory in which LaTeX finds them.

You may have to move the file with the extension .bst into a directory
in which BibTeX finds them.

Run the file aipcheck.tex through LaTeX. This will test your LaTeX
installation for completeness and reports any problems noticed.


Producing layout samples:
=========================

The file sample.tex can be used to produce samples for each proceeding
layout. For this it will interactively ask two questions when you run
through LaTeX:

 - any class options?    (you will normally answer with <RETURN> here
			  unless you wish to test some of the options
			  explained in aipguide.tex )

 - which layout style?   (here you have to name the layout style you
			  wish to produce, for example "6x9". A list
			  of choices will be presented at this point.)


So the way to proceed is

 latex sample              % answer questions
 bibtex sample		   % make bibliography
 latex sample		   % for fixing references
 latex sample		   % fix more references (depending on the
			   % style this might be necessary)

view or print sample.dvi



Precompiled samples (might not be distributed):
-----------------------------------------------

The distribution might contain a subdirectory or an archive with
various samples produced with the above method. 

sample-6s.ps   Postscript generated from sample.tex showing the AIP
	       6x9 proceedings layout

sample-8s.ps   Postscript generated from sample.tex showing the AIP
	       8.5x11 proceedings single column layout

sample-8d.ps   Postscript generated from sample.tex showing the AIP
	       8.5x11 proceedings double column layout




Further Notes:
==============

To run correctly the graphics package needs to be properly customized
for the site. This should be the case on standard installations; if
not, review cfgguide.tex and the graphics package documentation.

The natbib.sty is a powerful citation system which offers several
commands beside \cite. It comes with its own documentation, e.g.,
natnotes.tex.


Support for the ARLO journal:
=============================

ARLO has withdrawn LaTeX support for their journal, but the code is available
and can be extracted if the corresponding line in the installation file
(aipparms.ins and) is uncommented. There is no support for that journal.
