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
    style = <NDEx UUID of network to use for styling networks created>


The NDEx UUID needed for **style** can be obtained by uploading the :code:`style.cx` file found under
the :code:`data/` directory of this repository. NOTE: The network needs to be uploaded to the same
server as defined in **style** :code:`public.ndexbio.org` is NDEx_ production. Also the network needs
to be visible to the **user**

**Example configuration file**

.. code-block::

    [ncipid_dev]

    user = joe123
    password = somepassword123
    server = dev.ndexbio.org
    style = 86f63bf8-1b48-11e9-a05d-525400c25d22


Needed files
------------

SIF files can be found here: ftp://ftp.ndexbio.org/NCI_PID_EXTENDED_BINARY_SIF_2016-06-09-PC2v8%20API/

Files for :code:`--genesymbol, --networkattrib, --loadplan` flags can be found here:

https://github.com/coleslaw481/ndexncipidloader/tree/master/data


Usage
-----

For information invoke :code:`loadndexncipidloader.py -h`

**Example usage**

This example assumes a valid configuration file and needed files for :code:`--genesymbol, --networkattrib, --loadplan` are in the
current working directory and the SIF files are in the :code:`sif/` directory

.. code-block::

   loadncipidloader.py --genesymbol gene_symbol_mapping.json --loadplan loadplan.json --networkattrib netattrib.tsv <sif dir>


Via Docker
~~~~~~~~~~~~~~~~~~~~~~

**Example usage**

This example assumes files for :code:`--genesymbol, --networkattrib, --loadplan` are in the
current working directory, SIF files are in the :code:`sif/` directory, and a configuration
file has been created in current working directory and named :code:`conf`

.. code-block::

   docker run -v `pwd`:`pwd` -w `pwd` coleslawndex/ndexncipidloader:0.1.0 loadndexncipidloader.py --conf conf --genesymbol gene_symbol_mapping.json --loadplan loadplan.json --networkattrib netattrib.tsv sif

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _NDEx: http://www.ndexbio.org
