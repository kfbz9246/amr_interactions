This repo provides the code for the Preventative Veterinary Medicine paper 'Inferring the interaction structure of resistance to antimicrobials' - doi.org/10.1016/j.prevetmed.2018.02.007

A docker container with the environment required to run this code can be found at: kfbz9246/analysis_environment:10.1016-j.prevetmed.2018.02.007

The 'amr_db' directory contains the code for building the sqlite database that the analysis code queries.
	- 'create_narms_database.py' is the program that does the database creation
	- the 'database_interface' directory is a package implementing an sqlalchemy interface to the sqlite database that the analysis code requires.
	
The 'data_analysis' directory contains the code for all the analyses in the paper.
	- 'interactions_analysis_generalization.py' calculates and outputs interaction odds ratios and p-values. 
