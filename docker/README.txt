

NDEx NCI-PID Content Loader Docker Image
========================================

This image contains an installation of the Python package ndexncipidloader
from: https://github.com/coleslaw481/ndexncipidloader


To run:

   loadndexncipidloader.py

For help add -h:

   loadndexncipidloader.py -h

Files needed for --genesymbol, --loadplan, and --networkattrib can
be found in the /ndexncipidloader directory in this Docker image or
by visiting: https://github.com/coleslaw481/ndexncipidloader/data.

Also in the locations above is a style.cx which is an NDEx CX network that
contains styling for NCI-PID networks that can be uploaded to NDEx to get a
UUID to use for styling (its a required field in ~/.ndexutils.conf file)





