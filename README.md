PACTS -- Program for the Analytical Calculation of Transmission Spectra
---

This program uses the analytical formulae derived in Heng & Kitzmann (2017; https://arxiv.org/abs/1702.02051) 
in order to compute model transmission spectra. Line lists are obtained from ExoMol (http://exomol.com/).

Cross-sections
--------------
The cross sections are obtained from exo-mol (http://exomol.com/). In particular, they are assumed to go 
from 5000 cm^-1 to 29999 cm^-1 (i.e., go from 330 nm to 2000 nm), and have resolution of 0.1 cm^-1. Note 
that you have to download the cross-sections for your molecules from the webpage, and put them in the 
`cross_section` folder as `MOLECULE_TEMPERATURE.sigma` files (e.g., `H2O_1200.sigma` for the 
cross-section of water for a T = 1200 K atmosphere).
