=======
History
=======

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
