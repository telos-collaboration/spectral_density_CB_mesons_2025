# Analysis code for (arXiv:####.#####)

This is the analysis code built on top of **lsdensities** code (
<a href="https://github.com/LupoA/lsdensities"> GitHub repository </a>) to
reproduce the results in **([arXiv:####.#####][paper])**.

## Authors

Niccolò Forzano, Ho Hsiao, Fabian Zierler, Ed Bennett, Luigi Del Debbio, Ryan C. Hill,
Deog Ki Hong, Jong-Wan Lee, C.-J. David Lin, Biagio Lucini, Alessandro Lupo,
Maurizio Piai, Davide Vadacchino.

## Cloning the code

* Clone the repo including all submodules:
      git clone --recurse-submodules https://github.com/nickforce989/spectral_density_CB_mesons_2025.git


## Set up environments

* Download this code
* From the data release at https://doi.org/##.####/zenodo.########
  * Download ``chimera_data_reduced.hdf5``, and place it in ``input_correlators/``

* To automate find gevp plateaus, a x86-64 version of Snakemake must be used. To set this up, create a new x86-64 Conda environment with Snakemake, using
      
      conda create -n snakemake_x86 -c conda-forge -c bioconda snakemake
  
  and activate it
 
      conda activate snakemake_x86
 
  and then install Mamba
  
      conda install -c conda-forge 'mamba<2.0.0'


* Then, create the conda environment in terminal with conda installed:

      conda env create -f environment.yml
  
  with the caveat that if you're using an Apple silicon CPU then you need to use Conda 24.3 or later, and specify ```--platform osx-64```
  in your ```conda env create``` call.


* Once the environment is created, you can active the it:

      conda activate analysis-env4

* Then, install ``julia`` using conda

      conda install -c conda-forge julia

* This code uses a submodule: ```CB_autocorrelation_decay_constant```.
  Please, follow the instruction before using it, build the environment:
  
      conda env create --name wall_decay_constant --file=CB_autocorrelation_decay_constant/environment.yml



## Code usage

* The whole analysis can be done automatically:

   * To reproduce all the plots and results in the tables present in the paper, please run
     ``reproduce_everything.sh``. The results will be found in ``plots/`` and ``tables/``.
   * If ``latex`` is detected in the system, the full workflow will be ran. Otherwise, CSVs files will be produced in ``CSVs/`` and JSONs in ``JSONs/``. The tables will be generated in ``tables/``. 
        * In such a case, the output CSV and JSONs files can be used in a later moment to produce plots. This can be
        achieved by running ``run_plots.sh`` .

   * Make sure that the HDF5 file containing all the data is present in the  directory ``input_correlators/``.


* The analysis takes ~24 hours on the ```Sunbird``` supercomputer, where the CPU specifics are:
  ```2x Intel(R) Xeon(R) Gold 6148 CPU @ 2.40GHz with 20 cores each```.

## Acknoledgement

The flow_analysis code in ```topologies/flow_analysis``` has been based on the following <a href="https://github.com/edbennett/flow_analysis/"> GitHub repository </a>.

## License

[GPL](https://choosealicense.com/licenses/gpl-3.0/)


[paper]: https://arxiv.org/abs/####.#####
