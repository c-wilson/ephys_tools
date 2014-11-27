Electrophysiology Tools
====

This will be a comprehensive toolkit to combine data from electrophysiology recordings and Voyeur. Currently it is built
to work with raw data files recorded from spikeGL, but it can be easily updated to work with the future systems (think 
OpenEphys).

It is designed with the following design objectives in mind:
  
1.  <b>Data at rest should be transparent, easy to understand, and as minimally processed as is possible.</b> This allows 
users to parse data the way that makes sense for analysis at analysis run-time.
2.  Data-at-rest is stored using a powerful and open format once processed (HDF5).
3.  No substructure is assumed, and continuously recorded data are stored continuously (with the same time base) at rest.
Information to recreate substructures (ie trials) are stored with the data.
4.  Only two files should be required once data are processed: a raw data file and a processed data file. All metadata 
should be included within these files and propagated through the processing system.  
  

It is designed with the following packages:

ephys_tools.data_handling
----

This package contains packages for:  
  
1.  <u>ephys_tools.data_handling.pre-sort</u> parses data for sorting with the KlustaSuite of sofware (citation).  
2.  <u>ephys_tools.data_handling.post-sort</u> parses data into a structure that is useful for analysis. (See below for 
file structure).  
3.  <u>ephys_tools.data_handling.data_classes</u> module provides OO python interface for loading and manipulating data 
created by the above packages.  
4.  <u>ephys_tools.data_handling.exporters</u> package providing modules to export HDF5 file data to other file formats 
(ie Matlab). Not currently implemented and may not be necessary.  
  
<b>This also contains a file_structure_readme.md as a reference to the file structure.</b>


ephys_tools.analysis
----
Analysis modules. TBD.
