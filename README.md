# Pu-seq.R version 1.0
## Description
From the 4 bin-count files (count.csv: Pol-delta & Pol-episilon of Watson & Crick strands), wig format files containing track data of Pol-delta and Pol-epsilon usage are outputted. The bedgraph format file containing data of replication origins can be generated by setting option ?ori. The outputted wig files are the bedgraph files can be inputted into [ Integrative Genomics Viewer (IGV)](https://www.broadinstitute.org/igv/ "I").
## Comand
R < Pu-seq.R --vanilla --slave --args --prefix [name] --df  pol-delta.watson.count.csv --dr  pol-delta.crick.count.csv  --ef pol-epsilon.watson.count.csv --er pol-epsilon.watson.count.csv -ma1 [int] -ma2 [int] (-ori)
## Options
-prefix: The name of the output files; For example, if WT�h is entered here, the name of outputted file will be like �gWT_Pol-delta.crick.wig�h (required)  
-df, -dr, -er, -ef: The location and file names of the inputting csv files  (required)   
-ma1: the �gmoving average for ratios (default: 3)  
-ma2: the �gmoving average for origins, only meaningful when -?ori is set (default: same as --ma1)  
 -ori: activating origin detection and generation of the bedgraph file containing data on origin positions and efficiencies  
*Do not move -?args option it needs to be before all the options above.*
## Example
R < Pu-seq.R --vanilla --slave --args --prefix strain_name --df  pol-d.e1.f-w300.count.csv --dr  pol-d.e1.r-w300.count.csv  --ef pol-e.e1.f-w300.count.csv --er pol-e.e1.r-w300.count.csv -ma1 3 -ori