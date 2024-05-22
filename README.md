ChemDampPersson
===

The python script to compute the adsorbate induced contribution to line-width parameters $\Delta A$, and surface plasmon width 	$\gamma$, for chemical interface damping according to the method described by Persson, B. N. J., [Polarizability of small spherical metal particles: influence of the matrix environment](https://www.sciencedirect.com/science/article/abs/pii/003960289390865H). *Surface Science* **1993**, *281*, 153-162.

The latest version of the code is available at [https://github.com/blevine37/ChemDampPersson](https://github.com/blevine37/ChemDampPersson).

 ## Citation
 Publications arising from the use of this code are asked to please cite the forthcoming paper:

 B. Ostovar, et al.  The Role of the Plasmon in Interfacial Charge Transfer. (2024).
 
 ## Requirements
The script requres `Python3`, `numpy`, `scipy` and `matplotlib`. 

## Input Parameters
Enter the following parameters in the header section of `Calc_ChemDampPersson.py` file:  
* Planck's constant (eV.s)
* Fermi energy (eV) of core metal
* Energy (eV) of adsorbate induced resonance or virtual state (See [Persson Fig. 2](https://www.sciencedirect.com/science/article/abs/pii/003960289390865H))
* Observed position (eV) of Surface Plasmon Resonance (SPR)
* $\Gamma$ line width (eV) of adsorbate induced resonance or virtual state 
* Distance (Å) between dynamic image plane of and center of mass of orbital 
* Fermi velocity of core metal (Å/s)
* Q for p<sub>x</sub> or p<sub>y</sub> orbital of adsorbate (See [Persson Page 156](https://www.sciencedirect.com/science/article/abs/pii/003960289390865H))
* Carrier concentration of core metal (Å<sup>-3</sup>)
* Number of adsorbate per unit surface area (Å<sup>-2</sup>)
* Radius (Å) of particle
* Bulk dielectric constant of adsorbate matrix

## Output
The script gives the normal and tangential contributions of the adsorbate to the line-width parameters $\Delta A$, and surface plasmon width $\gamma$. The script also generates the plots of the function *J*(ω), and function Im $\alpha$(ω). For details see [Personn's paper](https://www.sciencedirect.com/science/article/abs/pii/003960289390865H).      

 ## Test
The test directory contains the applications of script to Ag@CO matrix, and reproduces the results given in [Persson's paper](https://www.sciencedirect.com/science/article/abs/pii/003960289390865H).

 ## Further Readings
For further details, please see the following publication.

1. Persson, B. N. J., [Polarizability of small spherical metal particles: influence of the matrix environment](https://www.sciencedirect.com/science/article/abs/pii/003960289390865H). *Surface Science* **1993**, *281*, 153-162.

## Contact
For any question, contact Ben Levine (ben.levine@stonybrook.edu) or A. Mehmood: (arshad.mehmood@stonybrook.edu).

## License
See the LICENSE file.
