Pre-sort package
====

This contains both scripts to process data for sorting with KlustaSuite software and templates for parameter (and probe)
files that are used to run these scripts. The main routine is process.py, which will transfer binary files (.bin) files into structured documents (raw.kwd). It will
also combine voyeur data into these documents.

Currently data are processed in the following hierarchy:  
><b>Sessions</b>  
Sessions are made up of multiple recs, each with their own raw.kwd file. 
>> <b>Recs</b>  
Recs may consist of one or many runs that are concatenated into a single raw.kwd file. 
<u>Recs are clustered as a single entity.</u>
>>> <b>Runs</b>  
Runs consist of individual bin files from the acquisition system and voyeur HDF5 files.  

When complete, a typical session will be made up of multiple rec files, but all runs will be included in the recs.

klustasuite
----
The primary dependancy of this preprocessing routine is the klustasuite. This software To install as of 7 Dec 2014 on linux:  

1. Download and install anaconda python distribution  
2. Once installed, you must create your python environment. The commands below will create this environment, named 'klusta'
<pre><code>$ conda create -n klusta python=2.7 --yes  
$ conda install -n klusta scipy pandas=0.12 pytables=3.0 pyqt setuptools pip cython nose ipython-notebook matplotlib --yes  
$ conda install -n klusta numpy=1.8 --yes</code>  

3. You must activate the klusta environment for the process.py script to run successfully. To activate or deactivate the 
environment:
<pre><code>$ source activate klusta</code></pre>

4. follow the installation instructions at: https://github.com/klusta-team/example




