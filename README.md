# RNALig
Scoring function to predict binding affinity of RNA-drug complexes. 
## Instructions to use: 
First save the PDB files in your working directory. Then install all the packages and then run the Features_extraction oython script. During running it will give two options: 1- Want to extract features for all the PDBs
         2- Want to extract features for a specific PDB.
User can select option according to the rquirement. The script will save the results in the form of .csv file in working directory.
Then use Random Forest Regressor python script to predict the affinity. Need to give path of the working directory or working folder.
It will give .csv file with predicted binding affinities by addining one column in the same feature file.
To check the contribution of features in binding affinity prediction , in Random Forest Regressor the code of equation is also attached. You have to run it saperatly and you will get the equation showing the intercept of RFR and feature contributions.
 
