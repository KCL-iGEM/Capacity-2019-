# Capacity
Synthetic biology relies heavily upon the integration of mathematics and engineering principles to design new biological systems with novel functionalities. Alongside our wet lab studies, we have created two software tools which serve as an informative resource for users to not only investigate the size limitations of current viral vectors but also equip them with a tool to investigate how capsid architecture may be manipulated to accommodate for a candidate therapeutic gene. 

CapsidOptimizer

Calculate new icosahedral geometries to design novel theoretical viral capsids with optimised packaging capacities for the delivery of a select gene associated with a rare genetic disease. 



CapsidBuilder 

Investigate the feasibility of novel viral capsid construction by evaluating the protein expression levels from wet lab constructs using translation efficiency measurements.  These measurements will be used to select one of our constructs that may be to build novel theoretical viral capsids. 


Development of our software tools was achieved using Python 3.7. Data for capsid dimensions was retrieved from the VIPERdb database and interactive web-based subplots were designed using the plotly Python library. The results are output to the user as a combination of detailed tables printed within the shell of the user’s IDE, and also through web-based subplots demonstrating highlighted components of data.

To run our programs, all files need to be downloaded, P1.py contains variable lists, the functions for our calculation are found in P2.py Run the CapsidOptimizer and CapsidBuilder, from their respective py. files. The plotly Python library must also be installed to display web-based interactive subplots of the data. Once the files are run, user input for a one of 57 candiate genes is required.

