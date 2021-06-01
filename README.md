# Read and compensate BD Influx FCS operator

##### Description

`read_and_compensate_fcs` operator transforms and compensates FCS files from the
BD Influx to Tercen datasets.

##### Usage

Input projection|.
---|---
`documentId`        | is the documentId (document can be a single FCS file, or a zipped set of FCS files
                      and a compensation matrix)

Output relations|.
---|---
`filename`          | character, the name of the FCS file
`channels`          | numeric, one variable per channel in the FCS file

##### Details

The operator transforms FCS files from the BD Influx into Tercen dataset and 
directly transforms the channels according to the compensation matrix. If the 
document is a ZIP file containing a set of FCS files, the operator extracts the 
FCS files and tranforms them into Tercen datasets. 

For compensation, if a `.csv` file is included in the .ZIP file, it will apply this
as compensation matrix. Otherwise it will try to use the compensation matrix in the
header of the .FCS file. 

The Flow Cytometry Standard is a data file standard for the reading and writing 
of data from flow cytometry experiments.
