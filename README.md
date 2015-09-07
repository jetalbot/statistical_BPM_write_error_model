# Bit Patterned Media Write Model#

A statistical write error model for exchange coupled composite and single layer bit patterned media (BPM) islands. 

##What is it?#

It is simple, but realistic model of island switching that captures the essential features of data storage. The model enables realistic servo requirements to be established for BPM for a given bit error rate (BER). It avoids the assumptions of simpler models without resorting to micromagnetic simulations over huge populations of islands. The model can use a full head field distribution, which allows for head field assymmetry and non-linearity, to calculate write errors. Distributions of island position, geometric and magnetic parameters can also be included for a given distrubution function. 

The basic routines for calculating demagnetising factors, energy barriers and coercivity have been separated and are contained in the /Basic-Routines/ folder. These form the building blocks of the entire model, but can be used for simpler calculations. Examples of how to use these functions in various combinations are contained in the /Basic-Examples/ folder .   


##Further information#

The [wiki](https://github.com/jetalbot/statistical_BPM_write_error_model/wiki) gives details on how the functions of the model work, how to use the function seperately and how to combine the functions.
