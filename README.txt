This folder presents the MATLAB scripts created to performed the topology synthesis, parametric study, and image processing

%% Topology Synthesis
The codes to perform the topology synthesis consist of:

• graphGrammar.m: 	  	Main script to call the topology synthesis. Used to reproduce the graph grammar algorithm as created by Zimmerman.
• matlabPTU source code:  	The source code of matlabPTU by Walker. This code has been adjusted and appended in order to be used in the topology synthesis. 
			  	The code consists of multiple class definitions. The main class as provided is origami.m, which is used to create origami topologies, crease patterns, and to calculate its folded states.
• Filters: 		  	Functions used to filter the resulting graphs for isomorphism and symmetry

%% Parametric Study

%main functions:
• paramatric_study_main.m:	Main script of the parametric study. This script calls the scripts to generate the valid geometries, generate rigidly foldable paths, and evaluate these paths separately to ensure good memory management
• getValidGeometries.m:		Script that generates all valid geometries inside the given parameter range
• generateRFPaths.m		Generates all possible rigidly foldable paths
• evaluatePaths.m		Evaluates all generated rigidly foldable paths for feasibility
• Analysis.m			Evaluates all feasible paths based on the proposed criteria in order to find the optimal solution for the given design case

%% Other scripts
• imageProcessing.m		Used to determine the end effector coordinates of the physical prototype in the experiment
• maxAuxDimensions.m		Used to calculate the maximum auxiliary dimensions

