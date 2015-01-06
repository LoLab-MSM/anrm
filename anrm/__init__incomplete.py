"""
ANRM: Apoptosis Necroptosis Reaction Model
========================================

ANRM is a family of models of the TNF apoptosis and necroptosis pathways written
using the Python modeling framework PySB. Model variants focus on exploring
alternative hypotheses for the involvement of Bid in suppressing RIP1 mediated
necroptosis.

Documentation is available in the docstrings, doc/ directory and online at
http://earm.readthedocs.org .

Modules
-------
::

 albeck_modules   --- components for albeck_* models
 lopez_modules    --- components for lopez_* models
 shen_modules     --- components for chen_*, cui_* and howells models
 shared           --- shared constants and macros for all models

 everything else (including mito.*)
                  --- the models

Scripts
-------
::

 estimate.py      --- simple parameter estimation using simulated annealing
 model_specs.py   --- display number of rules/odes/params for all models
 test_models.py   --- minimal test to ensure models contain no blatant errors

"""
