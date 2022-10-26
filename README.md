# Simulation of optogenetic perturbation of membrane surface charge  

This repository hosts the code for the article: _“Spatiotemporal dynamics of membrane surface charge regulates cell polarity and migration”_. DOI:[10.1038/s41556-022-00997-7](https://doi.org/10.1038/s41556-022-00997-7). The programs collectively model different optogenetic perturbations of inner membrane surface charge inside an excitable network topology. please see the _Methods_ section of the paper for details of the model. 

Please read the [citation](#restrictions) details below if you want to use/incorporate/modify a part of this repository in your research. 

## File Details

- Fig7BC.m generates the nullclines for  Figure 7B,C

- For the rest, the following naming convention is used:

        - (figure_name).m is the main MATLAB function that generates the figure specified by the figure name. 
        - (figure_name)_sdefile.m is the supporting function file. It contains the respective SDE model.

    For example, Fig7DE.m generates the figures  Figure 7D,E.

- The auxiliary functions:

    - SDE_euler.m
    - SDE_split_sdeinput.m


## Requirements

The MATLAB 2019a (or later versions) should be installed in the system to run these programs. [SDE Toolbox](http://sdetoolbox.sourceforge.net/) was used to develop these programs. 


## Restrictions

This program is a free software (please see [License](#license) for details). However, if you are using this code in your work, please cite our work as:


> **Tatsat Banerjee, Debojyoti Biswas, Dhiman Sankar Pal, Yuchuan Miao, Pablo A. Iglesias, Peter N. Devreotes** _“Spatiotemporal dynamics of membrane surface charge regulates cell polarity and migration”_, Nature Cell Biology, 24, 1499-1515 (2022). DOI:[10.1038/s41556-022-00997-7](https://doi.org/10.1038/s41556-022-00997-7).


## Authors

The code was primarily developed by Debojyoti Biswas and Pablo A. Iglesias (Department of Electrical and Computer Engineering, Johns Hopkins University, Baltimore, MD, USA). 


## License 

Copyright © 2022 Tatsat Banerjee, Debojyoti Biswas, Dhiman Sankar Pal, Yuchuan Miao, Pablo A. Iglesias, and Peter N. Devreotes. 

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 
