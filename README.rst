===========================
NDEx NCI-PID content loader
===========================


.. image:: https://img.shields.io/pypi/v/ndexncipidloader.svg
        :target: https://pypi.python.org/pypi/ndexncipidloader

.. image:: https://img.shields.io/travis/coleslaw481/ndexncipidloader.svg
        :target: https://travis-ci.org/coleslaw481/ndexncipidloader

.. image:: https://readthedocs.org/projects/ndexncipidloader/badge/?version=latest
        :target: https://ndexncipidloader.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status


Python application that loads NCI-PID data into NDEx_

This tool parses


Tools
-----

* **loadndexncipidloader.py** -- Loads NCI-PID networks from SIF files to NDEx_

Dependencies
------------

* `ndex2 <https://pypi.org/project/ndex2>`_
* `ndexutil <https://pypi.org/project/ndexutil>`_
* `biothings_client <https://pypi.org/project/biothings-client>`_
* `requests <https://pypi.org/project/requests>`_
* `pandas <https://pypi.org/project/pandas>`_


Compatibility
-------------

* Python 3.3+

Installation
------------

.. code-block::

   git clone https://github.com/coleslaw481/ndexncipidloader
   cd ndexncipidloader
   make dist
   pip install dist/ndexncipidloader*whl


Configuration
-------------

The **loadndexncipidloader.py** requires a configuration file in the following format be created.
The default path for this configuration is :code:`~/.ndexutils.conf` but can be overridden with
:code:`--conf` flag.

**Format of configuration file**

.. code-block::

    [<value in --profile (default ndexncipidloader)>]

    user = <NDEx username>
    password = <NDEx password>
    server = <NDEx server(omit http) ie public.ndexbio.org>


**Example configuration file**

.. code-block::

    [ncipid_dev]

    user = joe123
    password = somepassword123
    server = dev.ndexbio.org


Needed files
------------

Several steps of processing are needed to generate the SIF files.

**TODO:** Add a script this tool to generate these SIF files

1) Download **owl.gz** files from: ftp://ftp.ndexbio.org/NCI_PID_BIOPAX_2016-06-08-PC2v8-API/

2) Gunzip **.owl** files

3) Download paxtools.jar (http://www.biopax.org/Paxtools/) (requires Java 8+)

4) Run the following command for each **.owl** file putting **.sif** files in a directory by themselves:

.. code-block::

    java -jar paxtools.jar toSIF <INPUT OWL FILE> <INPUT OWL FILE WITH .owl replaced with .sif> \
             "seqDb=hgnc,uniprot,refseq,ncbi,entrez,ensembl" \
             "chemDb=chebi,pubchem" -useNameIfNoId -extended


Files for :code:`--genesymbol, --networkattrib, --style, --loadplan` flags can be found under the **data/** subdirectory
of this repository.

* **netattrib.tsv** appears to be a tab delimited export of `this excel file <https://github.com/NCIP/pathway-interaction-database/blob/master/download/NCI-Pathway-Info.xlsx>`_
* **gene_symbol_mapping.json** appears to have come from a previous run of ncipid processing with `this script <https://github.com/ndexbio/ndexutils/blob/master/ndexutil/ebs/ebs2cx.py>`_
* **style.cx** Original style to use
* **normalizedstyle.cx** Newer normalized style

Usage
-----

For information invoke :code:`loadndexncipidloader.py -h`

**Example usage**

This example assumes a valid configuration file and needed files for :code:`--genesymbol, --networkattrib, --loadplan` are in the
current working directory and the SIF files are in the :code:`sif/` directory

.. code-block::

   loadncipidloader.py --genesymbol gene_symbol_mapping.json --style style.cx --loadplan loadplan.json --networkattrib netattrib.tsv <sif dir>


Via Docker
~~~~~~~~~~~~~~~~~~~~~~

**Example usage**

This example assumes files for :code:`--genesymbol, --networkattrib, --loadplan` are in the
current working directory, SIF files are in the :code:`sif/` directory, and a configuration
file has been created in current working directory and named :code:`conf`

.. code-block::

   docker run -v `pwd`:`pwd` -w `pwd` coleslawndex/ndexncipidloader:0.1.0 loadndexncipidloader.py --conf conf --genesymbol gene_symbol_mapping.json --style style.cx --loadplan loadplan.json --networkattrib netattrib.tsv sif


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _NDEx: http://www.ndexbio.org
