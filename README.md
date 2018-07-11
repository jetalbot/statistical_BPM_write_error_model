# Bit Patterned Media Write Model#

A statistical write error model for exchange coupled composite and single layer bit patterned media (BPM) islands. 

## What is it?

It is simple, but realistic model of island switching that captures the essential features of data storage. The model enables realistic servo requirements to be established for BPM for a given bit error rate (BER). It avoids the assumptions of simpler models without resorting to micromagnetic simulations over huge populations of islands. The model can use a full head field distribution, which allows for head field assymmetry and non-linearity, to calculate write errors. Distributions of island position, geometric and magnetic parameters can also be included for a given distrubution function. 

The basic routines for calculating demagnetising factors, energy barriers and coercivity have been separated and are contained in the /Basic-Routines/ folder. These form the building blocks of the entire model, but can be used for simpler calculations. Examples of how to use these functions in various combinations are contained in the /Basic-Examples/ folder .   


## Further information

The [wiki](https://github.com/jetalbot/statistical_BPM_write_error_model/wiki) gives details on how the functions of the model work, how to use the functions seperately and how to combine the functions.

## Usage

The statistical BPM write model and its associated components are free to use. We ask that any research published using this model cites one of the following papers:

[1] J. Kalezhi, B. D Belle, and J. J Miles. Dependence of Write-Window on Write Error Rates in Bit Patterned Media. Magnetics, IEEE Transactions on, 46(10):3752–3759, 2010.

[2] J. Kalezhi, J.J. Miles, and B.D. Belle. Dependence of Switching Fields on Island Shape in Bit Patterned Media. Magnetics, IEEE Transactions on, 45(10):3531– 3534, 2009.

[3] J. Kalezhi, S. J. Greaves, Y. Kanai, M. E. Schabes, M. Grobis, and J. J. Miles. A statistical model of write-errors in bit patterned media. Journal of Applied Physics, 111(5):053926, 2012.

[4] J. E. Talbot, J. Kalezhi, C. Barton, G. Heldt and J. Miles, "Write Errors in Bit-Patterned Media: The Importance of Parameter Distribution Tails," in IEEE Transactions on Magnetics, vol. 50, no. 8, pp. 1-7, Aug. 2014.

[5] J. E. Talbot, J. Kalezhi and J. Miles, "Determining the Anisotropy of Bit-Patterned Media for Optimal Performance," in IEEE Transactions on Magnetics, vol. 51, no. 11, pp. 1-4, Nov. 2015.

*
*
*
*





