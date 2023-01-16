# STARlight
This is the new repository for the STARlight Monte Carlo, which simulates ultra-peripheral collisions of relativistic heavy ions.   It simulates two-photon 
and coherent and incoherent photoproduction of vector mesons.  It also includes an interface to DPMJet to be able to simulate other, general photoproductio reactions.

The code was described in [STARlight: A Monte Carlo simulation program for ultra-peripheral collisions of relativistic ions](arXiv:1607.03838) by Spencer R. Klein, Joakim Nystrand, Janet Seger,Yuri Gorbunov and Joey Butterworth, Computer Physics Communications, 212, 258 (2017) also available at arXiv:1607.03838.

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

10. <a id="extra"></a> __Extra on setting up__ `slight.in`: There are a few test files in the `config` folder. These test files if they will be used must be renamed to `slight.in` and moved or copied to the installation directory: `build`. These files can also be edited to suite the simulation conditions as desired by the user, useful comments on the different parameters can be found in these test files. More detailed information on the parameters can be found in the [PDF documentation](Readme.pdf)


