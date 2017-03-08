## Convert the format of the delta file 
 
input : deltas repository containing one file per plate 
output: deltas fits file readable by lya_auto_correlation 
subtract(optional): substract the mean of the deltas and correct the factor due to the substraction of the continuum (eq 6 and 7, Baustista (2017))

comand: 
convert_deltas.py -d deltas -o delta.fits --substract 

## Computation of the Lya auto correlation function 

command :  lya_auto_correlation -f delta.fits -o xiA.data  

