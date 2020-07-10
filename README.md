# Quartet Sampling #

![alt text](https://github.com/FePhyFoFum/quartetsampling/blob/master/doc/logo.png)

Quartet Sampling is a method to analyze molecular phylogenies by calculating branch support using repeated sampling of quartets.  Quartet Sampling differs from other support methods by combining a set of of tests into a single, efficient framework to address phylogenetic discordance.  QS is particularly useful for very large and data-sparse alignments or large phylogenomic datasets.

## Latest Update ##
Version 1.3.1 - Adds RESULTS.node.count.csv that provides counts and representative topologies for the three quartet arragements.  This is useful when QD is low in determining the more common discordant option.

## Usage ##

For complete details, see the manual for Quartet Sampling (https://github.com/FePhyFoFum/quartetsampling/blob/master/quartetsampling.pdf).

## NEW! Visualization ##

Shuiyin Liu has created a new R-script that can visualize QS results:

[https://github.com/ShuiyinLIU/QS_visualization](https://github.com/ShuiyinLIU/QS_visualization)

Please reference his URL separately, if you use this visualization script.

## How to Cite ##

James B Pease, Joseph W Brown, Joseph F Walker, Cody E Hinchliff, Stephen A Smith. 2018. Quartet Sampling distinguishes lack of support from conflicting support in the green plant tree of life. American Journal of Botany. 105(3): 385–403. doi:10.1002/ajb2.1016

Link to publication: http://dx.doi.org/10.1002/ajb2.1016

Please also include the URL (https://www.github.com/fephyfofum/quartetsampling).

**Be sure also to cite RAxML-ng, RAxML, PAUP, or IQ-TREE** (depending on which engine you use).

## Contributors ##

* James B. Pease - [www.peaselab.org](http://www.peaselab.org)
* Stephen A. Smith [blackrim.org](http://blackrim.org)
* Cody Hinchliff 
* Joseph W. Brown
* Joseph F. Walker

## License ##

This file is part of 'quartetsampling'

'quartetsampling' is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

'quartetsampling' is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with Foobar.  If not, see (http://www.gnu.org/licenses/).
