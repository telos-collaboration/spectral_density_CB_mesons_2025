# Analysis code for (arXiv:####.#####)

This is the analysis code built on top of **lsdensities** code (
<a href="https://github.com/LupoA/lsdensities"> GitHub repository </a>) to
reproduce the results in **([arXiv:####.#####][paper])**.

## Authors

Niccolò Forzano,
Ho Hsiao,
Fabian Zierler,
Ed Bennett,
Luigi Del Debbio,
Ryan C. Hill,
Deog Ki Hong,
Jong-Wan Lee,
C.-J. David Lin,
Biagio Lucini,
Alessandro Lupo,
Maurizio Piai,
Davide Vadacchino.

## Requirements

- Conda, for example, installed from [Miniforge][miniforge]
- [Snakemake][snakemake], which may be installed using Conda
- LaTeX, for example, from [TeX Live][texlive]


## Setup

1. Install the dependencies above.
2. Clone this repository including submodules
   (or download its Zenodo release and `unzip` it)
   and `cd` into it:

   ```shellsession
   git clone --recurse-submodules https://github.com/telos-collaboration/spectral_density_CB_mesons_2025
   cd spectral_density_CB_mesons_2025
   ```

3. Download the data and metadata from the [data release][datarelease]:

   1. Download the file `chimera_data_reduced.h5`,
      and place it in the `input_correlators` directory
   2. Download the file `metadata.zip`,
      and unzip it at root level to create the `metadata` directory.
   3. Download the file `raw_data.zip`,
      and unzip it at root level to create the `raw_data` directory.

## Running the workflow

The workflow is run using Snakemake:

``` shellsession
snakemake --cores 1 --use-conda
```

where the number `1`
may be replaced by
the number of CPU cores you wish to allocate to the computation
(or `all` to use all available CPU cores).

Snakemake will automatically download and install
all required Python and Julia packages.
This requires an Internet connection.
If you are running in an HPC environment where you would need
to run the workflow without Internet access,
these may be prepared by
running the following commands on a login node with Internet access:

``` shellsession
snakemake --cores 1 --sdm conda --conda-create-envs-only
snakemake --cores 1 --use-conda cb_julia_instantiated
```

Once this is complete,
the remainder of the workflow will run without Internet access.

The full workflow takes around 24 hours to run on
a single node of the ```Sunbird``` supercomputer,
which has 2x Intel(R) Xeon(R) Gold 6148 CPUs @ 2.40GHz with 20 cores each.

## Output

Output plots, tables, and definitions
are placed in the `assets/plots` and `assets/tables` directories.

Output data assets are placed into the `data_assets` directory.

Intermediary data are placed in a number of other subdirectories,
including:

- `intermediary_data`
- `lsd_out`
- `JSONs`
- `CSVs`


## Reusability

This workflow is relatively tailored to the data
which it was originally written to analyse.
Extending the analysis by adding more ensembles has not been tested;
in addition to adding the relevant raw data and metadata,
changes to the workflow definition files
(under `workflow/`)
and to the underlying Python code
(under `src/` and `lsd_out/`),
would be necessary.

## Acknowledgment

Portions of this workflow are derived from
[the workflow][antisymmetric-workflow]
for the publication
[Lattice investigations of the chimera baryon spectrum in the Sp(4) gauge theory][antisymmetric-paper].

This workflow builds on other work by members of this collaboration and others,
much of which is referenced by Git submodule or by environment specification.
This includes:

- [flow_analysis](https://github.com/edbennett/flow_analysis)
- [format_multiple_errors](https://github.com/edbennett/format_multiple_errors)
- [MadrasSokal](https://github.com/fzierler/MadrasSokal)
- [Snakemake][snakemake]

## License

This workflow is shared under the [GNU General Public License 3.0][gpl].

[antisymmetric-paper]: https://doi.org/10.48550/arXiv.2311.14663
[antisymmetric-workflow]: https://zenodo.org/records/10929539
[datarelease]: https://doi.org/10.5281/zenodo.TODO_ZENODO_ID
[gpl]: https://choosealicense.com/licenses/gpl-3.0/
[miniforge]: https://github.com/conda-forge/miniforge
[paper]: https://doi.org/10.48550/arXiv.TODO_ARXIV_ID
[snakemake]: https://snakemake.github.io
[snakemake-conda]: https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html
[texlive]: https://tug.org/texlive/
