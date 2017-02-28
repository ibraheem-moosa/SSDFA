# SSDFA  
This is the repository for Scatter Search for DNA Fragment Assembly tool (SS-DFA).

# Usage
**java ScatterSearch FragmentFile ResultFile n m nHCIter thresholdWeight numOfGens diversityMeasure numOfRuns nSSIter**

Here numOfGens is the maximum number of generations executed by the Scatter Search algorithm, numOfRuns is the number of time Scatter Search algorithm is executed and the best result of all these runs is given as output and nSSIter is the maximum number of generations allowed by Scatter Search without any improvements. 
To disable the effect of the last parameter set it greater than or equal to numOfGens.

The output file will contain data from all the runs. Each output line of each run  contains:
RunId fitness totalOverlap ContigCount WallClockTime

# configuration
Default configuration is already hard-coded in the source code. The default (tuned) values are specified in the paper. Certain configurations can be overridden using a config.txt file. The file should be placed in the current directory (i.e. the directory from where the ScatterSearch class gets loaded).

The parameters "Selection strategy for fit individuals" and "Fitness function" can be changed by modifying the code. We made changes in the code directly when we experimented with different selection strategies and fitness functions. We do not recommend users changing these. However, if you must, you can go into the source code and make necessary changes.

To use traditional fitness function instead of our proposed one:
- Open up the code for Utility.fitness().
- Uncomment the line: //fitness = totalOverlap;
- Comment out the line right above it: fitness=((double)totalOverlap)/contigCount;

To use tournament selection instead of truncation selection for fit individuals:
- Open up the code for ScatterSearch.getFitIndividuals()
- Uncomment the line: //return getTournamentSelectedIndividuals(population, fitnessArray, n);
- Comment out the line right above it: return getTopRankedIndividuals(population, fitnessArray, n);

To disable exploitative mutation in hill-climbing (HC):
- Open up the code for LocalSearch.HillClimbing()
- Comment out the line: //int[] r = Utility.mutation(best, threshold, ov);
- Uncomment the line right above it: int[] r = Utility.mutation(best);

To run GAGE validation scripts you will have to compile and install MUMmer. **Go into the MUMmer directory and run make**. 
For this command to be successful there cannot be any spaces in the names of any directory in the path of the MUMmer directory. Also you may need to **install csh and some other dependencies**. 

Sample GAGE validation script run.
Go to gage-validation directory and run.
**sh getCorrectnessStats.sh ../Original/ACIN02000001.1.fasta ../FinalAssembly/en/ACIN02000001_Assembly_en_best.fa**

