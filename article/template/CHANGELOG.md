# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a
Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

## [v3.13c]
### Added
- Style for _ACS Bio. Med. Chem. Au_ (`abmcb8`)
- Style for _ACS Eng. Au_ (`aeacb3`)
- Style for _ACS Env. Au_ (`aeacc4`)
- Style for _ACS Mater. Au_ (`amacgu`)
- Style for _ACS Meas. Au_ (`amachv`)
- Style for _ACS Nanosci. Au_ (`anaccx`)
- Style for _ACS Org. Inorg. Au_ (`aoiab5`)
- Style for _ACS Phys. Chem. Au_ (`apcach`)
- Style for _ACS Polym. Au_ (`appccd`)
- Style for _JACS Au_ (`jaaucr`)

### Changes
- Minor adjustments from ACS

## [v3.13b]
### Added
- Style for _Acc. Mater. Res._ (`amrcda`)
- Style for _ACS Agr. Sci. Tech._ (`aastgj`)
- Style for _ACS Food Sci. Tech._ (`afsthl`)

## [v3.13a]
### Fixed
- Swap TOC height and width to give landscape output

## [v3.13]
### Added
- Style for _ACT ES&T Eng._ (`aeecco`)
- Style for _ACT ES&T Water_ (`aewcaa`)
- Changelog as separate file

### Changed
- Split `.ins` file out from `.dtx`

### Fixed
- Corrected standard TOC height and width (see #30 and #34)

## [v3.12a] - 2019-02-14
### Added
- Style for _ACS Mater. Lett._ (`amlcef`)

### Changed
- Print article titles in _J. Am. Chem. Soc._ bibliography style
- Use ISO date formats throughout

## [v3.12] - 2018-09-15
### Added
- Style for _ACS Appl. Electron. Mater._ (`aaembp`)
- Style for _ACS Appl. Polym. Mater._ (`aapmcd`)

### Fixed
- Error if explicit `and others` is used in `.bib` file (issue #29)

## [v3.11b] - 2018-07-12
### Changed
- Maximum number of authors printed in _ACS Nano_

### Fixed
- Print abstract text at full size
- Spacing when `and others` is used in `.bib` file (issue #29)

## [v3.11a] - 2018-02-05
### Changed
- Style changes for Supporting Information

## [v3.11] - 2018-01-10
### Added
- Style for _ACS Appl. Energy Mater._ (`aaemcq`)
- Style for _ACS Appl. Nano Mater._ (`aanmf6`)

### Changed
- Use `chemformula` in demo file
- Remove redundant comments about e-TeX
- Include article titles in all manuscript bibliographies

## [v3.10i] - 2017-05-18
### Changed
- Updates for the _Biochemistry_ style
- Updates for the _J. Org. Chem._ style
- Updates for _Org. Lett._ style

### Fixed
- Issue printing DOI in bibliographuy

## [v3.10h] - 2017-01-21
### Added
- Style for _ACS Earth Space Chem._ (`aesccq`)
- Style for _ACS Sensors_ (`aidcbc`)
- Style for _ACS Infect. Dis._ (`ascefj`)

### Changed
- Updates for the _J. Chem. Theory Comput._ style
- Updated text used for Supplementary Information

## [v3.10g] - 2017-01-14
### Added
- Style for _ACS Biomater. Sci._ (`abseba`)

### Fixed
- Issue in style for _Inorg. Chem._

## [v3.10f] - 2016-09-07
### Changed
- Updated text for Supplementary Information
- Updates for the _Macromolecules_ style

## [v3.10e] - 2016-08-31
### Changed
- Updates for the _Macromolecules_ style

## [v3.10d] - 2016-06-17
### Changed
- Improved DOI support

## [v3.10c] - 2016-05-10
### Added
- Style for _ACS Omega_ (`acsodf`)
- Style for _ACS Energy Lett._ (`aelccp`)

### Changed
- Added `achemso-` to all style configuration file names

## [v3.10b] - 2016-01-20
### Changed
- Updated text for Supplementary Information
- Split demonstration files from main source

### Fixed
- Define `\acs@tocentry@text` globally (issue #18)

## [v3.10a]
### Changed
- Updates for the _Chem. Mater._ style
- Updates for the _Chem. Rev._ style

## [v3.10] - 2015-04-07
### Changed
- Print affiliation information in the title block

### Fixed
- Include `\section*` and `\subsection*` to TOC

## [v3.9b] - 2015-03-23
### Changed
- Add T1 encoding to demo file (issue #14)
- Superscript references for _J. Agric. Food Chem._ style (issue #16)

### Fixed
- Allow `super=true` to work after loading `natbib` (issue #15)

## [v3.9a] - 2015-03-12
### Changed
- Print DOIs in all cases in _ACS Central Sci._ style

## [v3.9] - 2015-02-01
### Added
- Option `doi`
- Version information to BibTeX log

## [v3.8n] - 2015-01-18
### Added
- Style for _ACS Central Sci._ (`acscii`)

### Changed
- Style changes for _ACS Nano_
- Remove keywords from _J. Phys. Chem._ styles
- TOC entry size for _J. Phys. Chem. Lett._

## [v3.8l] - 2014-08-23
### Changed
- Allow for case where fax number is given with no phone number

## [v3.8k] - 2014-08-18
### Fixed
- Format date correctly if given for `in press` BibTeX entries

## [v3.8j] - 2014-06-12
### Changed
- Update article title and keyword requirements for all journals

### Fixed
- Configuration file name for _Biochemistry_ (`bichaw`)

### Removed
- Style for_Biotechnol. Prog._

## [v3.8i] - 2014-05-14
### Changed
- Style update for _Chem. Res. Toxicol._
- Include keywords for _ACS Appl. Mater. Interfaces_ style

## [v3.8h] - 2014-03-31
### Changed
- Style updates for _Chem. Res. Toxicol._

## [v3.8g] - 2014-03-06
### Changed
- Formatting for _J. Phys. Chem._ styles

### Fixed
- Ensure that `mciteplus` patch is safe even if package loading is bypassed

## [v3.8f] - 2014-01-23
### Changed
- Update following `mciteplus` change

## [v3.8e] - 2014-01-08
### Changed
- Change TOC entry size for _J. Med. Chem._
- Include titles for article in _ACS Appl. Mater. Interfaces_ style

## [v3.8d] - 2013-10-04
### Fixed
- Typo in keyword printing

## [v3.8c] - 2013-09-20
### Changed
- Updated _ACS Photonics_ style

## [v3.8b] - 2013-09-15
### Changed
- Alter TOC printing slightly
- Update demo to reflect TOC changes

### Fixed
- Print abstract/TOC/keywords in correct order for _J. Phys. Chem. Lett._

## [v3.8a] - 2013-09-08
### Changed
- Update bibliography style for _J. Phys. Chem._ derivatives

### Fixed
- Minor style fixes

## [v3.8] - 2013-08-23
### Added
- Command `\latin`
- Style for _ACS Applied Materials & Interfaces_
- Style for _ACS Photon._ (`apchd5`)
- Style for _ACS Sustainable Chem. Eng._ (`ascecg`)
- Style for _ACS Synth. Biol._ (`asbcd6`)
- Style for _Environ. Sci. Technol. Lett._ (`estlcu`)

### Changed
- Print keywords after abstract not after title
- Drop loading any font files
- Extend author symbols to 99 authors
- Print keyword and abbreviations in sections

### Fixed
- Behaviour of `\citenum`
- Minor style fixes

## [v3.7h] - 2013-06-07
### Fixed
- Remove extraneous `*` in format.doi function (issue #13)

## [v3.7g] - 2013-04-13
### Changed
- Update manuscript types for _ACS Macro Lett._
- Print article title for _J. Phys. Chem._ papers

### Fixed
- Set TOC size for _ACS Macro Lett._ (issue #10)
- Add DocStrip guards for _ACS Macro Lett._ (issue #11)

## [v3.7f] - 2013-03-07
### Changed
- Load `mciteplus` after `natbib`

### Fixed
- `etalmode` option

## [v3.7e] - 2013-02-13
### Changed
- Drop use of `mathptmx`

## [v3.7d] - 2012-09-12
### Fixed
- Typo in the `.bst`

## [v3.7c] - 2012-08-30
### Changed
- Update for _Ind. Eng. Chem. Res._ style

## [v3.7b] - 2012-07-25
### Changed
- Update _J. Phys. Chem. Lett._ style to truncate author list

### Fixed
- Step version (issue #8)

## [v3.7a] - 2012-05-08
### Changed
- Print article titles in _Langmuir_ style

### Fixed
- Relationship between `usetitle` and `articletitle`

## [v3.7] - 2012-04-30
### Added
- Option `chaptertitle`
- Style for _ACS Med. Chem. Lett._ (`amclct`)

### Changed
- Rename `usetitle` to `articletitle` (to match `biblatex` `chem-acs` style)

### Fixed
- No author truncation for _J. Chem. Theory Comput._ (issue #7)

## [v3.6] - 2012-04-30
### Changed
- Drop `cleveref`

## [v3.5k] - 2012-04-30
* 8aca9c0 Various incredibly niggly things for J. Nat. Prod.

## [v3.5j] - 
* 9a5cb7f Flexible size for TOC entry box

## [v3.5i] - 
* 1649618 Print keywords for J. Phys. Chem. submissions

## [v3.5h] - 
* 93811a3 Give schemes, etx. 'S' numbers
* e8d62c8 Use 'S<num>' for SI refs.
* 7478673 A slight improvement for SI
* b34463a Formatting tweak for J. Chem. and Eng. Data.

## [v3.5g] - 
* ee474b2 Another .bst typo

## [v3.5f] - 
* 920dc5c One-token typo in .bst affecting 'and others' entries

## [v3.5e] - 
* d44219d Checking new formats: minor adjustment
* 5834bfb Settings for new journals (hopefully)
* 83c9401 New ACS journals for 2011:  - New docstrip guards Remove J. Comb. Chem.
* 1be7600 Remove engine line (not really needed)

## [v3.5d] - 
|\  
| * 6955482 Fix bug in make.ordinal
* | 2431abf Added tag v3.5c for changeset beb648a09793
* |   1da12a1 Resolve merge issues
|\ \  
| |/  
|/|   
| * 9feff46 Resolve update issues
| ## [v3.5a] - 
| * 73bbe01 Make \mciteSubRef work properly (fixes #6) Update Org. Lett. formatting (fixes #5) Bump version number
* | 3513189 Added tag v3.5b for changeset 8ee4387af241
* | bd472bb Bug fix for nalefd.cfg
* | c84d4bf Added tag v3.5b for changeset 0b1acfd46340
* | fae782e Added tag v3.5a for changeset 20b1bf86cfe9
* | 166a701 Print abstract for JACS Communications (fixes #3)
* |   6ff6c54 Merge changes
|\ \  
| |/  
| * b44e931 Add keywords to Nano Lett. output
| ## [v3.5] - 
| * 9df91e0 Note behaviour of booktitles/chapters
| * 2c926dd Lots of improvements following the style guide
* | 4e174ec Move biblabel stuff to end of class/package (fixes #4)
|/  
* 4267fa3 Corrected formatting based on Style Guide: patents and techreports
* 4cfd1f1 Use method for defining \endmcitethebibliography that is not caught out by LaTeX2e's \newcommand restriction on "\end..." functions
* 1dae139 Step version number (for tomorrow)
* 60b25c5 A couple of fixes from the style guide≈y
*   a67ff26 Clean up head (what is up?)
|\  
| * 55e1aae Adjust position of that new test reference
* | e7a81e7 Logic reversal in BibTeX stlyes => do the same in package/class code
* | b90db69 Move position of year in Biochemistry-style articles
* | db54abb Completely revised bst files: rewritten from scratch Still need to check versus ACS Style Guide
* | 7c33384 Added 'patent' record to demo database
|/  
* 97f4038 Add "inbook" reference to test doc (for testing updates to BibTeX style)
* cb1920d Use natbib mechanism to alter \@biblabel

## [v3.4g] - 
* 4815f00 Added ACS Chem. Neurosci. style

## [v3.4f] - 
* 15cba28 Only define bibnote counter if not set up elsewhere Style file for J. Chem. Ed.

## [v3.4e] - 
* 78dfe94 Biomacromolecule style improvements again

## [v3.4d] - 
* a48eb33 Improved BibTeX style for reports

## [v3.4c] - 
* ad71645 More improvements for Biomacromolecules

## [v3.4b] - 
* 015b3fc Add hyperref option so things work properly with cleveref

## [v3.4a] - 
* e0bf02e Better style file for Biomacromolecules

## [v3.4] - 
* c7e9034 New advice on installation Switch from varioref to cleveref for \ref support New TOC printing routine Include fixes from new version of cite package in natmove

## [v3.3g] - 
* a4871fe Advice on overriding class options added

## [v3.3f] - 
* 68d0914 Style fix for Acc. Chem. Res.

## [v3.3e] - 
* b221e15 Minor style fixes

## [v3.3d] - 
* 8638068 Bug fix in meta data routines

## [v3.3c] - 
* 1d50b41 Add advice on mismatching class and BibTeX style files REmove spurious comma when two author names given Add J. Phys. Chem. Lett. style

## [v3.3b] - 
* d0f63c7 Bug fix for affiliation symbol printing

## [v3.3a] - 
* 874a7b3 New options usetitle and etalmode for controlling bibliography

## [v3.3] - 
* 457f0cb New layout option for two column printing Better logic for \altaffiliation to reduce duplicate printing Improve keyval code (use less xkeyval functions)

## [v3.2f] - 
* 99425ef Better installation advice Move bibliography set up to correct bug with biochem style Improved footnote symbol macro

## [v3.2e] - 
* e787692 Alteration to the method for printing the abstract, to allow delayed print

## [v3.2d] - 
* 72b0817 Fix addition of star to lead author

## [v3.2c] - 
* 537aa6f Fix bug with email option New \fax and \phone meta-data options Use optional argument for \title as running header suggestion

## [v3.2b] - 
* 14963fd Added email option to set printing of e-mail address on first page

## [v3.2a] - 
* b813cf6 Bug fix for duplicate institute test

## [v3.2] - 
* ddd8fe8 New functions to alter formatting (sections, section numbers, abstract) Environment to format TOC entry added Some internal refactors

## [v3.1a] - 
* b000d9f Bug fix for duplicate affiliation detection

## [v3.1] - 
* f998e90 Produce print-like layout for JACS Communications Add \alsoaffiliation as a complement to \altaffiliation Convert \suppinfo and \acknowledgement to environments Package loading reordered Meta data system rearranged

## [v3.0a] - 
* 5de99e0 Skip printing footnotes when only one institution is needed

## [v3.0]
### Changed
- Second re-write, converting to a class and giving much tighter integration
  with ACS submission system

## [v2.0]
### Changed
- Re-write of package by Joseph Wright

## [v1.0]
### Added
- Initial release of package by Mats Dahlgren

[Unreleased]: https://github.com/josephwright/achemso/compare/v3.13bc..HEAD
[v3.13b]: https://github.com/josephwright/achemso/compare/v3.13b...v3.13c
[v3.13b]: https://github.com/josephwright/achemso/compare/v3.13a...v3.13b
[v3.13a]: https://github.com/josephwright/achemso/compare/v3.13...v3.13a
[v3.13]: https://github.com/josephwright/achemso/compare/v3.12a...v3.13
[v3.12a]: https://github.com/josephwright/achemso/compare/v3.12...v3.12a
[v3.12]: https://github.com/josephwright/achemso/compare/v3.11b...v3.12
[v3.11b]: https://github.com/josephwright/achemso/compare/v3.11a...v3.11b
[v3.11a]: https://github.com/josephwright/achemso/compare/v3.11...v3.11a
[v3.11]: https://github.com/josephwright/achemso/compare/v3.10i...v3.11
[v3.10i]: https://github.com/josephwright/achemso/compare/v3.10h...v3.10i
[v3.10h]: https://github.com/josephwright/achemso/compare/v3.10g...v3.10h
[v3.10g]: https://github.com/josephwright/achemso/compare/v3.10f...v3.10g
[v3.10f]: https://github.com/josephwright/achemso/compare/v3.10e...v3.10f
[v3.10e]: https://github.com/josephwright/achemso/compare/v3.10d...v3.10e
[v3.10d]: https://github.com/josephwright/achemso/compare/v3.10c...v3.10d
[v3.10d]: https://github.com/josephwright/achemso/compare/v3.10b...v3.10c
[v3.10b]: https://github.com/josephwright/achemso/compare/v3.10a...v3.10b
[v3.10a]: https://github.com/josephwright/achemso/compare/v3.10...v3.10a
[v3.10]: https://github.com/josephwright/achemso/compare/v3.9b...v3.10
[v3.9b]: https://github.com/josephwright/achemso/compare/v3.9a...v3.9b
[v3.9a]: https://github.com/josephwright/achemso/compare/v3.9...v3.9a
[v3.9]: https://github.com/josephwright/achemso/compare/v3.8n...v3.9
[v3.8n]: https://github.com/josephwright/achemso/compare/v3.8m...v3.8n
[v3.8m]: https://github.com/josephwright/achemso/compare/v3.8l...v3.8m
[v3.8l]: https://github.com/josephwright/achemso/compare/v3.8k...v3.8l
[v3.8k]: https://github.com/josephwright/achemso/compare/v3.8j...v3.8k
[v3.8j]: https://github.com/josephwright/achemso/compare/v3.8i...v3.8j
[v3.8i]: https://github.com/josephwright/achemso/compare/v3.8h...v3.8i
[v3.8h]: https://github.com/josephwright/achemso/compare/v3.8g...v3.8h
[v3.8g]: https://github.com/josephwright/achemso/compare/v3.8f...v3.8g
[v3.8f]: https://github.com/josephwright/achemso/compare/v3.8e...v3.8f
[v3.8e]: https://github.com/josephwright/achemso/compare/v3.8d...v3.8e
[v3.8d]: https://github.com/josephwright/achemso/compare/v3.8c...v3.8d
[v3.8c]: https://github.com/josephwright/achemso/compare/v3.8b...v3.8c
[v3.8b]: https://github.com/josephwright/achemso/compare/v3.8a...v3.8b
[v3.8a]: https://github.com/josephwright/achemso/compare/v3.8...v3.8a
[v3.8]: https://github.com/josephwright/achemso/compare/v3.7h...v3.8
[v3.7h]: https://github.com/josephwright/achemso/compare/v3.7g...v3.7h
[v3.7g]: https://github.com/josephwright/achemso/compare/v3.7f...v3.7g
[v3.7f]: https://github.com/josephwright/achemso/compare/v3.7e...v3.7f
[v3.7e]: https://github.com/josephwright/achemso/compare/v3.7d...v3.7e
[v3.7d]: https://github.com/josephwright/achemso/compare/v3.7c...v3.7d
[v3.7c]: https://github.com/josephwright/achemso/compare/v3.7b...v3.7c
[v3.7b]: https://github.com/josephwright/achemso/compare/v3.7a...v3.7b
[v3.7a]: https://github.com/josephwright/achemso/compare/v3.7...v3.7a
[v3.7]: https://github.com/josephwright/achemso/compare/v3.6...v3.7
[v3.6]: https://github.com/josephwright/achemso/compare/v3.5k...v3.6
[v3.5k]: https://github.com/josephwright/achemso/compare/v3.5j...v3.5k
[v3.5j]: https://github.com/josephwright/achemso/compare/v3.5i...v3.5j
[v3.5i]: https://github.com/josephwright/achemso/compare/v3.5h...v3.5i
[v3.5h]: https://github.com/josephwright/achemso/compare/v3.5g...v3.5h
[v3.5g]: https://github.com/josephwright/achemso/compare/v3.5f...v3.5g
[v3.5f]: https://github.com/josephwright/achemso/compare/v3.5e...v3.5f
[v3.5e]: https://github.com/josephwright/achemso/compare/v3.5d...v3.5e
[v3.5d]: https://github.com/josephwright/achemso/compare/v3.5c...v3.5d
[v3.5c]: https://github.com/josephwright/achemso/compare/v3.5b...v3.5c
[v3.5b]: https://github.com/josephwright/achemso/compare/v3.5a...v3.5b
[v3.5a]: https://github.com/josephwright/achemso/compare/v3.5...v3.5a
[v3.5]: https://github.com/josephwright/achemso/compare/v3.4g...v3.5
[v3.4g]: https://github.com/josephwright/achemso/compare/v3.4f...v3.4g
[v3.4f]: https://github.com/josephwright/achemso/compare/v3.4e...v3.4f
[v3.4e]: https://github.com/josephwright/achemso/compare/v3.4d...v3.4e
[v3.4d]: https://github.com/josephwright/achemso/compare/v3.4c...v3.4d
[v3.4c]: https://github.com/josephwright/achemso/compare/v3.4b...v3.4c
[v3.4b]: https://github.com/josephwright/achemso/compare/v3.4a...v3.4b
[v3.4a]: https://github.com/josephwright/achemso/compare/v3.4...v3.4a
[v3.4]: https://github.com/josephwright/achemso/compare/v3.3g...v3.4
[v3.3g]: https://github.com/josephwright/achemso/compare/v3.3f...v3.3g
[v3.3f]: https://github.com/josephwright/achemso/compare/v3.3e...v3.3f
[v3.3e]: https://github.com/josephwright/achemso/compare/v3.3d...v3.3e
[v3.3d]: https://github.com/josephwright/achemso/compare/v3.3c...v3.3d
[v3.3c]: https://github.com/josephwright/achemso/compare/v3.3b...v3.3c
[v3.3b]: https://github.com/josephwright/achemso/compare/v3.3a...v3.3b
[v3.3a]: https://github.com/josephwright/achemso/compare/v3.3...v3.3a
[v3.3]: https://github.com/josephwright/achemso/compare/v3.2f...v3.3
[v3.2f]: https://github.com/josephwright/achemso/compare/v3.2e...v3.2f
[v3.2e]: https://github.com/josephwright/achemso/compare/v3.2d...v3.2e
[v3.2d]: https://github.com/josephwright/achemso/compare/v3.2c...v3.2d
[v3.2c]: https://github.com/josephwright/achemso/compare/v3.2b...v3.2c
[v3.2b]: https://github.com/josephwright/achemso/compare/v3.2a...v3.2b
[v3.2a]: https://github.com/josephwright/achemso/compare/v3.2...v3.2a
[v3.2]: https://github.com/josephwright/achemso/compare/v3.1a...v3.2
[v3.1a]: https://github.com/josephwright/achemso/compare/v3.1...v3.1a
[v3.1]: https://github.com/josephwright/achemso/compare/v3.0a...v3.1
[v3.0a]: https://github.com/josephwright/achemso/compare/v3.0...v3.0a
