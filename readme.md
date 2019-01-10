# Anomaly Screening Applying Pseudowells (ASAP)

_FORCE Hackathon, Stavanger 18-19 September 2018._

### Done

* merge pseudo-well blocking with notebook 1
 
### To do

* add rock physics modeling (Per)
* add seismic analysis (Lukas): ONGOING (separate repo)
* clean up code (Sandeep)
* check everything (ALL)


### Additional notes

In `ASAP_notebook_1` I have included a first section that does a bit of introduction to the well datasets, reading in all 3 wells from Per Avseth's QSI dataset, creating a few stratigraphic markers and plottings a summary of all the key logs.

However during the hackathon we used a version of well 15/5-5 which Anders had in matlab format and I have saved as csv, changing the log names. As I remember, this is the "correct" version with a more accurate porosity log and fluid-replaced logs. It would be good however to have a more complete log suite with also water saturation and the insitu elastic logs. 

For example, these are the logs contained in qsiwell2.csv (equivalent to 15/5-5):

    for i in range(w2.columns.size):
       print('Log {} = {}'.format(i,w2.columns[i]))

    Log 0 = VP
    Log 1 = VS
    Log 2 = RHO_OLD
    Log 3 = GR
    Log 4 = NPHI
    Log 5 = RHO
    Log 6 = SWE
    Log 7 = SWX
    Log 8 = VPVS
    Log 9 = IP
    Log 10 = IS
    Log 11 = VSH
    Log 12 = RHOm
    Log 13 = RHOf
    Log 14 = PHIE

I have included all the main functions in a separate file (`asap_library.py`). We have to agree as to keep it like this (and include perhaps larger functions used by Lukas?) or integrate all the functions within the Jupyter notebooks.

I have also modified extensively Anders' code (essentially his `generate_wells`), without messing it up (I think). There is a few things a bit obscure, like the relations used to modify the elastic logs on the basis of saturation and porosity variations. These comes from the rock physics relations established by Per, so we need to include his analysis too. A plot would be enough as a start, but then we need to explicitly code it in the notebooks (and maybe create a direct link with the pseudowells generation).
