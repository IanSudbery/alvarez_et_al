# Pipeline and code for Alvarez-Benayas, Katsarou and Trasanidis et al #

## Dependency installation ##
Dependencies are listed in environment.yml. This can be used to build a conda environment. We recommend using mamba to do this, as it is significantly faster than conda.

```
$ mamba env create [-p PATH_TO_STORE_ENV|-n NAME_FOR_ENV] -f environment.yml
```

Next you will need to install `cgat-flow`, which is not yet available on conda. This can be achieved by cloning the github repository and installing it into your environment. 

```
$ conda activate ENV_NAME_OR_PATH
$ git clone https://github.com/cgat-developers/cgat-flow.git
$ cd cgat-flow
$ python setup.py develop
```

## Configuring the pipeline ##

.. TODO

## Running the pipeline ##

.. TODO

