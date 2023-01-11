# nLayerModel

**GENERAL INFO**

This repository contains all code that I created during my master thesis.  
Used Versions:
 - Python 3.9.12
 - Spyder 5.3.1
 - Comsol 5.2


 **DIRECTORY OVERVIEW**
 
The key directory of my thesis is ```modules```. It contains all the math to create the *nLayerModel*. 
- ```model.py``` contains all top level functions to setup the nLayerModel
- ```submodel.py``` is used to create and solve the matrix equation for all current loadings
- ```layer.py``` stores all the lowest level math that is the solution forms of the vector potential, flux densities and field strength
- the subdirectory ```plot``` is called for visualization

The directory ```analytics``` contains the analytical solutions of a 2,3 and 4 region model.  
It also stores some basic formulas that are used from other scripts

The directory ```launchers``` stores the comparison of analytic results and the solution of the *nLayerModel*.  
It also contains the comparison of the *nLayerModel* to other studies.

In the ```data``` directory all the fixed values for the 7 MW generator investigaed in my thesis are defined.  
It also holds the code to export variables to COMSOL Multiphysics.

The directory ```design``` stores most of the research that was done with the *nLayerModel*.  
Many scripts here are kind of messy due to a lot of on the fly changes.  
Noteworthy is the file ```n7_model.py``` which contains many functions and classes to define and investigate 7 layer models.

Finally the directory ```thesis``` contains most scripts and data that are used in my thesis.  

**SETUP for running all scripts in Spyder** 

To run scripts that are not inside the projects root directory,  
select the root directory in ```Tools -> PYTHONPATH manager -> add path```.  
With that all modules load from any script in the project.
