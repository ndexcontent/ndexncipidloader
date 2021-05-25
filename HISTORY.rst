=======
History
=======

5.0.1 (2021-05-25)
-----------------------

* Switched default layout (can be overridden with `--layout` flag) to `force-directed`
  since `force-directed-cl` may not work on all machines.

5.0.0 (2021-05-20)
-----------------------

* Fixed duplicate node issue by removing nodes and edges from a network if a family node, contains
  the node in its `memberlist`. Any edges are shifted to the family node with duplicates
  merged where possible.

4.0.0 (2020-11-04)
-------------------

* New default behavior: **force-directed-cl** layout is now applied on
  networks via py4cytoscape library and a running instance of Cytoscape.
  Alternate Cytoscape layouts and the networkx "spring" layout can be
  run by setting appropriate value via the new **--layout** flag

3.1.1 (2020-10-16)
-------------------

* Removed NODE_LABEL_POSITION discrete mapping from style since it is
  not compatible with CX 2.0

3.1.0 (2019-09-11)
-------------------

* Added **--disableshowcase** flag that lets caller disable showcasing of **NEWLY** added networks which is enabled by default.

* Added **--indexlevel** flag that lets caller set type of indexing performed on **NEWLY** added networks. Default is full indexing (all).

3.0.0 (2019-08-02)
-------------------

* Renamed command line tool from **loadndexncipidloader.py** to **ndexloadncipid.py** to be more consistent with other loaders. Since this is a breaking change bumped to version 3.0.0

* Added **--visibility** flag which lets caller dictate whether newly added networks are set to PUBLIC (default) or PRIVATE

* Removed parameter **--disablcitededgemerge** since the changes in 2.0.0 causes this to no longer have any effect

* Set default for **--paxtools** flag to be **paxtools.jar** which assumes the tool is in current working directory

2.0.0 (2019-07-16)
------------------

* Spring layout adjusted by increasing iterations

* Code now removes all neighbor-of edges with NO data migration. controls-state-change-of
  edges are removed if more informative edges exist. Any orphaned nodes resulting from
  the removal of these edges are also removed

1.6.0 (2019-07-09)
------------------

* Added *__iconurl* network attribute to all networks

* Added **interactome** to *networkType** network attribute for 'NCI PID - Complete Interactions' network

1.5.1 (2019-07-09)
------------------

* Renamed network attribute *type* to *networkType* to adhere to normalization specification

1.5.0 (2019-06-28)
------------------

* Fixed style.cx by removing view aspects that was causing networks to not render properly in cytoscape

1.4.0 (2019-06-13)
------------------

* Network PathwayCommons.8.NCI_PID.BIOPAX is now renamed
  to 'NCI PID - Complete Interactions' with alternate description.

1.3.0 (2019-06-12)
------------------

* Improved description in style.cx file (JIRA ticket UD-362)

1.2.0 (2019-06-11)
------------------

* Code now adds a citation attribute to every edge even if there is no value
  in which case an empty list is set (JIRA ticket UD-360)

* Added type network attribute and set it to ['pathway'] following normalization
  guidelines

1.1.0 (2019-06-10)
------------------

* Adjusted network layout to be more compact by reducing number of iterations in
  spring layout algorithm as well as lowering the value of scale (JIRA ticket UD-360)

1.0.2 (2019-05-24)
------------------

* Removed view references from cyVisualProperties aspect of style.cx file cause it was causing issues with loading in cytoscape

* Set directed edge attribute type to boolean cause it was incorrectly defaulting to a string

1.0.1 (2019-05-18)
------------------

* Renamed incorrect attribute name prov:wasDerivedBy to prov:wasDerivedFrom
  to adhere to normalization document requirements
 
1.0.0 (2019-05-16)
------------------

* Massive refactoring and first release where code attempts to behave as defined in README.rst

0.1.1 (2019-02-15)
------------------

* Updated data/style.cx by renaming Protein to protein and SmallMolecule
  to smallmolecule to match the new normalization conventions


0.1.0 (2019-02-15)
------------------

* First release
