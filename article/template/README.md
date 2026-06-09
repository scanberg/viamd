achemso - Support for submissions to American Chemical Society journals
=======================================================================

The `achemso` bundle provides a LaTeX class file and BibTeX
style file in accordance with the requirements of the American
Chemical Society (ACS). The files can be used for any documents,
but have been carefully designed and tested to be suitable for
submission to ACS journals.

The bundle also includes the `natmove` package. This package is
loaded by `achemso`, and provides automatic moving of
superscript citations after punctuation.

Installation
------------

The package is supplied in `.dtx` format and as a pre-extracted
`.zip` file, `achemso.tds.zip`. The later is most convenient for
most users: simply unzip this in your local `texmf` directory.
If you want to unpack the `.dtx` yourself, running `tex
achemso.dtx` will extract the package whereas `latex
achemso.dtx` will extract it and also typeset the documentation.

Typesetting the documentation requires a number of packages in
addition to those needed to use the package. This is mainly
because of the number of demonstration items included in the
text. To compile the documentation without error, you will
need the packages:
 - `array`
 - `booktabs`
 - `hypdoc`
 - `listings`
 - `lmodern`
 - `mathpazo`
 - `microtype`