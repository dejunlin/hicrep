======
HiCRep
======

HiCRep is a python implementation of the R package with the same name. It provides a tool for assessing the reproducibility of Hi-C datasets based on the algorithm published by `Yang et al. <https://pubmed.ncbi.nlm.nih.gov/28855260/>`_
In brief, data is first smoothed to reduce the effects of noise and bias and then stratified by contact distance to adjust for the distance dependency of Hi-C datasets. These strata are compared between samples to calculate the stratum-adjusted correlation coefficient (SCC). The SCC is a value between -1 and 1 and can be interpreted the same way you would a standard correlation. 


This implementation takes a pair of Hi-C data sets in Cooler format (.cool for single binsize or .mcool formultiple binsizes) and computes the HiCRep SCC scores for each pair of chromosomes between the two data sets. It can be run either as a command line tool of through a python API. 

Installation
============
To install hicrep you will need python version 3.7.6 or greater. Additionally, hicrep has the following dependencies:

- `NumPy <https://numpy.org/>`_
- `SciPy <https://www.scipy.org/>`_
- `cooler <https://cooler.readthedocs.io/en/latest/datamodel.html>`_
- `pandas <https://pandas.pydata.org/>`_
- `h5py <https://www.h5py.org/>`_
- `statsmodels <https://www.statsmodels.org/stable/index.html>`_

We recommend using pip to install hicrep, which will install the package along with all of its dependencies. 

.. code-block:: bash

   $ pip install hicrep

Once hicrep is installed, you can interact with it either through the command line or as a python module.

CLI Guide
=========

hicrep can be run from the command line with:

.. program:: hicrep
.. code-block:: shell

    hicrep cool1 cool2 outfile [OPTIONS] 

.. rubric:: Arguments

.. option:: cool1

    The first Hi-C contact matrix in .cool or .mcool format [required]

.. option:: cool2

    The second Hi-C contact matrix in .cool or .mcool format [required]

.. option:: outfile

    A path to an output file to write the SCC scores for each chromosome [required]

.. rubric:: Options

.. option:: -h , --help

    See a list of expected arguments to hicrep

.. option:: --binsize <BINSIZE>

    The binsize you want to use from a .mcool file. Default is -1, meaning that the inputs are .cool files with a single precomputed binsize

.. option:: --h, <h>

    Used to set the size of the sliding 2D window used for smoothing the contact matrix. Size will be 1 + 2 * h bins. We recommend the following values for h based on the resolution of your data:

    +------------+------------+
    | resolution |     h      | 
    +============+============+
    |    10kb    |     20     | 
    +------------+------------+
    |    25kb    |     10     | 
    +------------+------------+
    |    40kb    |      5     |
    +------------+------------+
    |   100kb    |      3     | 
    +------------+------------+
    |   500kb    |   1 or 2   | 
    +------------+------------+
    |    1MB     |   0 or 1   |
    +------------+------------+


.. option:: --dBPMax <DBPMAX>

    Only consider contacts at most this number of bp away from the diagonal. Defaults to -1, meaning the entire contact matrix is used. In general, we recommend not using the entire contact matrix because usually the last few diagonals have very few valid data in it for computing Pearson's correlation

.. option:: --bDownSample

    When set, hicrep will down sample the input with more contact counts to the the same number of counts as the other input with less contact counts. If not set, the input matrices will be normalized by dividing the counts by their respective total number of contacts.


API Guide
=========
You can also use hicrep as a python module. First, use the util function `readMcool`

.. code-block:: python

    from hicrep.utils import readMcool

to read a pair of `.mcool` files and specify the bin size to compute SCC with:

.. code-block:: python

    fmcool1 = "mydata1.mcool"
    fmcool2 = "mydata2.mcool"
    binSize = 100000
    cool1, binSize1 = readMcool(fmcool1, binSize)
    cool2, binSize2 = readMcool(fmcool2, binSize)

or a pair of `.cool` files with built-in bin size:

.. code-block:: python

    fcool1 = "mydata1.cool"
    fcool2 = "mydata2.cool"
    cool1, binSize1 = readMcool(fmcool1, -1)
    cool2, binSize2 = readMcool(fmcool2, -1)
    # binSize1 and binSize2 will be set to the bin size built in the cool file
    binSize = binSize1

then define the parameters for computing HiCRep SCC:

.. code-block:: python
    
    from hicrep import hicrepSCC

    # smoothing window half-size
    h = 1

    # maximal genomic distance to include in the calculation
    dBPMax = 500000

    # whether to perform down-sampling or not 
    # if set True, it will bootstrap the data set # with larger contact counts to
    # the same number of contacts as in the other data set; otherwise, the contact 
    # matrices will be normalized by the respective total number of contacts
    bDownSample = False


    # compute the SCC score
    scc = hicrepSCC(cool1, cool2, h, dBPMax, bDownSample)

which will result in a list containing the SCC scores for each chromosome available in the data set.
