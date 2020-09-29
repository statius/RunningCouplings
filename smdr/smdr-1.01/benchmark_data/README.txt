Users should typically not modify the file ReferenceModel.dat in this 
directory! It is used in the compilation of SMDR. The ReferenceModel.dat
will over-write the one in the main directory after: 

  make clean
  make

To make this ReferenceModel.dat, the procedure we followed is:

  ./calc_fit -int -o NEW_ReferenceModel.dat -e 1.0e-13

and answer the prompts according to the latest PDF data. Then, edit 
NEW_ReferenceModel.dat to eliminate obvious roundoff error in inputs 
(for example Mt_pole as 172.4 rather than 172.39999999999) and 
eliminate trailing 0's. Then place the resuling file in this directory 
as ReferenceModel.dat.

Older versions of ReferenceModel.dat will also be stored here, with names

  ReferenceModel_20**.dat

where 20** is the year.
