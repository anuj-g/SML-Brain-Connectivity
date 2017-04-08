----------------------------------------------------------------------------------------------------
Brain Network Analysis For Gender Classification
----------------------------------------------------------------------------------------------------

This repository contains all the code that are needed to run the experiments for analyzing gender differences in human connectomes.

Pre-requisites:
1. MATLAB
2. python 
3. scikit-learn and its related pre-requisites:http://scikit-learn.org/stable/
4. classify: Contains the code for classification
5. bootstrap_generic: Contains the code for bootstrapping procedure
6. ccf: code to calculate ccf
7. connectivity: Connectivity heatmap code
8. csv_writes : Write CSV file which can be used as input to classifier
9. edge_betweenness: Code to calculate edge between-ness centrality
10. efficiency: Code to calculate local efficiency of each node
11. normalization_scripts: Code to normalize
12. participation_coeff: Code for participation coeff


The code is structured as follows:
1. BCT: Brain connectivity toolbox(also here): We use this toolkit to calculate network science measures:https://sites.google.com/site/bctnet/
2. heatmaps: Toolkit for drawing heatmaps in MATLAB: http://www.mathworks.com/matlabcentral/fileexchange/24253-customizable-heat-maps
3. datasets: All the data. The MATLAB code assumes that the data set directory has 2 sub-directories: 'males' and 'females' where the corresponding .mat and .csv
   files are present.
4. To convert the raw csv files to .mat file do csv_to_mat('datasets/<set1|set2>')
5. The normalization of the matrices is done by:
   normalizematrix('datasets/set2')


Any of the network measures can be calculated by running the corresponding script and passing the path to the normalized data set and relevant arguments as shown
below as an example:

Experiment Notes: The following commands assume that the directory have two sub directories, 'males' and 'females' where the corresponding .mat or .csv files are present.

CCF Based on mean
>> ccf_normalized_errorbars('datasets/set2/normalized')

CCF Based on median
>> ccf_normalized_median('datasets/set2/normalized')

If you need to run the classification code just type:
python brain_classify_set1.py
python brain_classify_set2.py


Vivek Kulkarni and Jagat Sastry P.
