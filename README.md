# Analysis code for (arXiv:2405.01388)

This is the analysis code built on top of **lsdensities** code (
<a href="https://github.com/LupoA/lsdensities"> GitHub repository </a>) to
reproduce the results in **([arXiv:2405.01388][paper])**.

## Authors

Niccolò Forzano, Ho Hsiao, Fabian Zierler, Ed Bennett, Luigi Del Debbio, Ryan C. Hill,
Deog Ki Hong, Jong-Wan Lee, C.-J. David Lin, Biagio Lucini, Alessandro Lupo,
Maurizio Piai, Davide Vadacchino.


## Set up environment

* Download this code
* From the data release at https://doi.org/10.5281/zenodo.11048346
  * Download ``chimera_data_full.hdf5``, and place it in ``input_correlators/``
  * Download ``input_topology.zip``,
    and extract its contents into    ``input_topology/``
  * Download ``input_fit/``,
    and extract its contentse into ``input_fit/``.


* Then, create the conda environment in terminal with conda installed:

      conda env create -f environment.yml
  
  with the caveat that if you're using an Apple silicon CPU then you need to use Conda 24.3 or later, and specify ```--platform osx-64```
  in your ```conda env create``` call.


* Once the environment is created, you can active the it:

      conda activate analysis-env

* Then, install ``julia`` using conda

      conda install -c conda-forge julia

## Code usage

* The whole analysis can be done automatically:

   * Make sure that ``input_fit/``, ``input_topology/`` and ``input_correlators/`` contain the relevant files from the data release,
   as discussed above.
   * To reproduce all the plots and results in the tables present in the paper, please run
     ``reproduce_everything.sh``. The results will be found in ``plots/`` and ``tables/``.
   * If ``latex`` is detected in the system, the full workflow will be ran. Otherwise, CSVs files will be produced in ``CSVs/``. 
        * In such a case, the output CSV files in ``CSVs/`` can be used in a later moment to produce plots and tex tables. This can be
        achieved by running ``plots_and_tables.sh`` .

* Otherwise, each individual step can be achieved separately:

   * Make sure that ``input_topology/`` is full. To reproduce all the plots that are shown in the paper, run 
     ``bash run_plots.sh``.  The results will be found in ``plots/``.

   * Make sure that the HDF5 file containing all the data is present in the 
     directory ``input_correlators/`` before running ``bash run_spectral_densities.sh``.

     The spectral densities can be found and fitted:
      * If ``input_fit/`` has been filled (using the corresponding directory in https://doi.org/10.5281/zenodo.11048346)
        the fitting procedure will be applied to pre-reconstructed spectral densities.
      * If ``input_fit/`` is empty, the code will reconstruct from scratch the spectral densities and then
        to fit them. This procedure may take quite a long time.
  
     Either way, to do so, run ``run_spectral_densities.sh``. 

   * Make sure that the HDF5 file containing all the data is present in the  directory ``input_correlators/``. The 
     command ``bash run_plateaus.sh`` reproduces the GEVP plateaus by using the  metadata used for the analysis in the paper.



## Acknoledgement

The flow_analysis code in ```topologies/flow_analysis``` has been based on the following <a href="https://github.com/edbennett/flow_analysis/"> GitHub repository </a>.

## License

[GPL](https://choosealicense.com/licenses/gpl-3.0/)


[paper]: https://arxiv.org/abs/2405.01388
