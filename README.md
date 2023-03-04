# STARlight
This is the new repository for the STARlight Monte Carlo, which simulates ultra-peripheral collisions of relativistic heavy ions.   It simulates two-photon 
and coherent and incoherent photoproduction of vector mesons.  It also includes an interface to DPMJet to be able to simulate other, general photoproductio reactions.

The code was described in [STARlight: A Monte Carlo simulation program for ultra-peripheral collisions of relativistic ions](https://doi.org/10.1016/j.cpc.2016.10.016) by Spencer R. Klein, Joakim Nystrand, Janet Seger,Yuri Gorbunov and Joey Butterworth, Computer Physics Communications, 212, 258 (2017) also available at [arXiv:1607.03838](https://arxiv.org/abs/1607.03838).

## Installation Instructions
The following instructions illustrates the procedure to install and run STARlight in a *nix based environment:

1. Download the code package from GitHub to the desired location. You can either download the latest release on GitHub or Clone the repository into your desired location - e.g. `Desktop` as shown below:
    ```
    cd ~/Desktop
    git clone https://github.com/STARlightsim/STARlight.git
    ``` 
2. Switch into the folder where STARlight has just downloaded - `Desktop/STARlight`
    ```
    cd ~/Desktop/STARlight
    ```
3. Create an installation directory and __switch to this directory__ - e.g. `build`
    ```
    mkdir build
    cd ~/Desktop/STARlight/build
    ```
4. Setup this installation directory with `cmake`. *Note: Ensure that you are in the installation directory, and that you have `cmake` installed.*
    ```
    cmake ~/Desktop/STARlight
    ```
5. Compile the code using either `make` or `gmake`
    ```
    (g)make
    ```
6. This compilation will produce an executable: `starlight` in the installation directory. Confirm this before you proceed.
7. Setup the desired running condition in the input file: `slight.in`. *Ensure that you have a `slight.in` file in the installation directory - `build` - before you run the executable*.
The easiest way to setup `slight.in` is to start with a test file and edit [as described here](#extra). But For your first simulation, you can use the __default__ `slight.in` as shown below:
    ```
    cp ~/Desktop/STARlight/config/slight.in ~/Desktop/STARlight/build
    ```
8. Run the simulation:
    ```
    ./starlight
    ```
9. The output will be found in the installation directory with the default name: `slight.out`. This name can change depending on the `basefilename` set in your `slight.in`. For example, if your `basefilename` is `slightvoppPb_lhc`, your output file will be named as `slightvoppPb_lhc.out` in the installation directory.

10. <a id="extra"></a> __Extra on setting up__ `slight.in`: There are a few test files in the `config` folder. If these test files will be used, they must be renamed to `slight.in` and moved or copied to the installation directory: `build`. These files can also be edited to suit the simulation conditions desired by the user. Useful comments on the role of different parameters can be found in these test files and more detailed information on the parameters can be found in the [PDF documentation](Readme.pdf)


## HEPMC3 OutPut
The Official source tarball for installation can be found [here](http://hepmc.web.cern.ch/hepmc/index.html)

After unpacking the source tarball, compile, build and install the HEPMC3 package into the ```desired_installation_path``` following the instructions [here](http://hepmc.web.cern.ch/hepmc/building.html). 

Please take note of the ```desired_installation_path``` where you installed HEPMC3, as you would need this location to link HEPMC3 to STARlight.

To compile STARlight with HEPMC3 output enabled use:
```
cmake /pathto/starlight -DENABLE_HEPMC3=ON -DHepMC3_DIR=/pathto/hepmc3/hepmc3-install
```

For example: if your ```desired_installation_path``` is ```~/Desktop/STARLIGHT/HepMC3-3.2.5/build```, and your present working directory is the STARlight build directory, then use:
```
cmake .. -DENABLE_HEPMC3=ON -DHepMC3_DIR=~/Desktop/STARLIGHT/HepMC3-3.2.5/build/
```