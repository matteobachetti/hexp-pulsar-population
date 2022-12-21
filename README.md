# hexp-pulsar-population

## Requirements
The user needs to have [SIXTE](https://github.com/thdauser/sixte), [SIMPUT](https://github.com/thdauser/simput), and [HEASOFT](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/download.html) installed and working.

Moreover, they need to have an up-to date version of response files and SIXTE configurations for HEX-P and any other instrument they want to simulate.

Finally, they need some sort of source catalog and, possibly, a FITS image with valid WCS containing the diffuse emission from the host galaxy or, when relevant, the interstellar medium, and its spectrum in a valid Xspec `.xcm` file.

## Usage

This repository contains two Python scripts and two example Bash scripts:

+ `sixte_population.py`: this produces a SIMPUT file containing many sources. Each of the sources is split into two:
   
   1. A thermal spectrum (`diskpbb` in Xspec), optionally with aperiodical variability (a complicated PDS with various Lorentzian components)
   2. A hard cutoff spectrum, pulsed. Pulsations have random pulsed fractions (from 0 to 100%), and the pulse profile is a Von Mises profile with $\kappa$ parameter ranging from 1 to 20.
   
   Launch with 
   
   ```
   $ python sixte_population.py catalog_file.fits output_file.simput diffuse_emission_image.img diffuse_emission_spectrum.xcm 9.6e-12
   ```
   where 9.6e-12 is the flux of the diffuse emission in this case.
   
+ `run_sim_hexp.bash` executes the simulation, using reasonable defaults. Launch with (e.g. a 100 ks observation)

  ```
  $ bash run_sim_hexp.bash output_file.simput 100000
  ```
  creates files named like `output_file_let_evt.fits`, `output_file_het_1_evt.fits`, `output_file_het_2_evt.fits`
  
+ `time_source_from_sixte.py` reads the simulation, the SIMPUT catalog, and calculates the H test of all the pulsars to test whether they are detected. Launch with, e.g., for a 5 arcsecond extraction region,
  ```
  $ python time_source_from_sixte.py output_file.simput output_file_het_1_evt.fits 5
    ```
  This produces, for each pulsar, a folded pulsar profile and an HDF5 file containing the data of the pulse profile (can be read with an astropy Table as `table=Table.read("file.hdf5")`)
