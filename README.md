# RNA-ligand
Scoring function to predict binding affinity of RNA-drug complexes. 
This scoring function is built to predict the binding affinity of RNA-Drug complexes.
## Instructions to use: 
First save the PDB files in your working directory. Then use all the three features code RNA specific features, Drug specific features and Complex specific features.These codes will extract all the features in a seperate files in .csv format. Them merge all the files of features make a single file and save it in your working directory.
Run the Random Forest Regressor with your data. 
It will give .csv file with predicted binding affinities.
To check the contribution of features in binding affinity prediction , in Random Forest Regressor the code of equation is also attached. You have to run it saperatly and you will get the equation showing the intercept of RFR and feature contributions.
 
