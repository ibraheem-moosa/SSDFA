/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author pritom
 */
import java.io.*;
import java.util.*;
import suffixtree.*;
import java.util.*;

public class ScatterSearch {

    public static void main(String args[]) {

        if(args.length != 10) {
            System.out.println("Usage: java ScatterSearch FragmentFile ResultFile n m nHCIter thresholdWeight numOfGens diversityMeasure numOfRuns nSSIter");
            System.exit(1);
        }
        String strResultFile = args[1];

        long fileReadStartTime = System.currentTimeMillis();
        
        GeneralizedSuffixTree tree = ReadAndStoreFragments(args[0]);
        int[][] overlapArray = tree.allPairSuffixPrefix();
        long fileReadEndTime = System.currentTimeMillis();

        System.out.println("Fragments, config read and store time: " + (fileReadEndTime - fileReadStartTime) + " ms");

        //
        // SIZE is total number of fragments in the tree, divided by 2. Because
        // Each fragment is kept twice in the tree; once in forward direction
        // another one as RC (reverse complement)
        //
        
        int SIZE = tree.getStringCount() / 2;

        /*
        int n = 30;
        int m = 20;
        int nHCIter = 500;
        double thresholdWeight = 0.20;
        int numOfGens = 200;
        String diversityMeasure = "PDistance";
        */

        int n = Integer.parseInt(args[2]);
        int m = Integer.parseInt(args[3]);
        boolean noDiversity = (m == 0);
        int popSize = n + m;
        int nHCIter = Integer.parseInt(args[4]);
        double thresholdWeight = Double.parseDouble(args[5]);
        int numOfGens = Integer.parseInt(args[6]);
        String diversityMeasure = args[7];
        int numOfRuns = Integer.parseInt(args[8]);
        int nSSIter = Integer.parseInt(args[9]);

        if (0 == diversityMeasure.compareTo("PDistance"))
            Utility.divMeasure = Utility.DivMeasure_PDistance;
        else
            Utility.divMeasure = Utility.DivMeasure_HamDistance;

        long totalOverlap = 0; 
        
        //
        // We set the threshold to be thresholdWeight fraction of mean fragment length
        //
        
        long totalLength = 0;
        for (int i = 0; i < SIZE; i++) {
            totalLength += tree.getString(i).length();
        }
        
        long threshold = (long)(thresholdWeight * totalLength / SIZE);
        System.out.println("Overlap threshold (Auto-tuned):  " + threshold);

        //
        // Setup numOfRuns random numbers as seeds for 10 runs
        //
        Utility.rng.setSeed(10);
        
        long[] runRandSeeds = new long[numOfRuns];
        for (int i = 0; i < runRandSeeds.length; i++)
            runRandSeeds[i] = Utility.rng.nextLong();

        double[] fitnessArray = new double[popSize * popSize];
        double[] diversityArray = new double[popSize * popSize];

        // best individual from all runs
        int[] BEST = null;
        double fitnessBest = Double.NEGATIVE_INFINITY;

        for (int runId = 1; runId <= numOfRuns; runId++)
        {
            long startTime = System.currentTimeMillis();

            System.out.println("RUN " + runId);
            System.out.println("-------------");

            // Set random number seed
            Utility.rng.setSeed(runRandSeeds[runId - 1]);

            int[][] seededPopulation = generateSeededPopulation(n, SIZE, overlapArray, nHCIter, threshold);
            System.out.println("Generated Seeded Population");

            int[][] diversePopulation = generateDiversePopulation(m, SIZE, seededPopulation);
            System.out.println("Generated Diverse Population");
            seededPopulation = null;
            int[][] population = new int[popSize * popSize][];
            for(int i = 0; i < popSize; i++)
                population[i] = diversePopulation[i];
            diversePopulation = null;
            
            // fitnessRunBest contains the best fitness of this run
            double fitnessRunBest = Double.NEGATIVE_INFINITY;
            // best individual of this run
            int[] runBEST = null;

            // Do hill climbing on initial population
            for (int i = 0; i < popSize; i++) {
                population[i] = LocalSearch.HillClimbing(overlapArray, population[i], nHCIter, threshold);
                fitnessArray[i] = Utility.fitness(population[i], overlapArray, threshold);
            }
             
            for (int i = 0; i < popSize; i++) {
                
                if (runBEST == null || fitnessArray[i] > fitnessRunBest) {
                    runBEST = population[i];
                    fitnessRunBest = fitnessArray[i];
                }
            }

            long runBESTTotalOverlap = Utility.overlapCount(runBEST, overlapArray, tree);
            long runBESTnContig = Utility.contigCount(runBEST, overlapArray, threshold);
            
            System.out.println("Initial Population Best Fitness: " + fitnessRunBest + " Total Overlap: " + runBESTTotalOverlap);
            
            int numOfGenWithNoUpgrade = 0;
            //////////////////////////////////////////////////////////////////////////////////////////////
            //MAIN LOOP STARTS FROM HERE
            /////////////////////////////////////////////////////////////////////////////////////////////
            for (int time = 1; time <= numOfGens; time++) {
                
                double fitnessGenBest = Double.NEGATIVE_INFINITY;
                
                // Extract diverse and fit individuals to generate parents
                // at time == 1 we already have it
                if(time != 1){
                    int[][] B = getFitIndividuals(population, fitnessArray, n);
                
                    //Evaluate Diversity, fitness has already been calculated
                
                    int[][] D = null; 
                    
                    if(!noDiversity){
                        for (int i = 0; i < population.length; i++) {
                            diversityArray[i] = Utility.diversity(population, population[i], population.length);
                        }
                        D = getDiverseIndividuals(population, diversityArray, m);
                    }

                    //Generating Parents for new Population
                    for (int i = 0; i < n; i++) 
                        population[i] = B[i];
                    for (int i = n; i < popSize; i++)
                        population[i] = D[i - n];
                    for (int i = popSize; i < population.length; i++)
                        population[i] = null;
                }
                //
                // For each pair, there will be 2 children generated. So, the size of newly generated
                // set of individuals would be: popSize * (popSize - 1)
                // Plus, the existing set size is popSize.
                //
                // So, the total size of new population: popSize * popSize
                //

                ArrayList<int[]> newIndividuals = new ArrayList<>(popSize * popSize);

                for(int i = 0; i < popSize; i++) {
                    for(int j = 0; j < i; j++) {
                        int[][] children = Utility.crossover(population[i], population[j]);
                        //
                        // Normal mutation operation
                        //
                        children[0] = Utility.mutation(children[0]);
                        children[1] = Utility.mutation(children[1]);
                        //
                        // Exploitative mutation operation
                        //
                        children[0] = Utility.mutation(children[0], threshold, overlapArray);
                        children[1] = Utility.mutation(children[1], threshold, overlapArray);
                        // Hill climbing
                        children[0] = LocalSearch.HillClimbing(overlapArray, children[0], nHCIter, threshold);
                        children[1] = LocalSearch.HillClimbing(overlapArray, children[1], nHCIter, threshold);

                        newIndividuals.add(children[0]);
                        newIndividuals.add(children[1]);
                    }
                }
                
                for(int i = 0; i < popSize; i++){
                    newIndividuals.add(population[i]);
                }
                
                newIndividuals.toArray(population);

                for (int i = 0; i < population.length - popSize; i++) {
                    fitnessArray[i] = Utility.fitness(population[i], overlapArray, threshold);
                    if (fitnessArray[i] > fitnessGenBest) {
                        fitnessGenBest = fitnessArray[i];
                    }
                    if (fitnessArray[i] > fitnessRunBest) {
                        runBEST = population[i];
                        fitnessRunBest = fitnessArray[i];
                    }
                }
                
                if(fitnessGenBest < fitnessRunBest){
                    numOfGenWithNoUpgrade++;
                }
                else{
                    numOfGenWithNoUpgrade = 0;
                    runBESTTotalOverlap = Utility.overlapCount(runBEST, overlapArray, tree);
                    runBESTnContig = Utility.contigCount(runBEST, overlapArray, threshold);
                }

                long endTime = System.currentTimeMillis();
                long duration = (endTime - startTime) +  (fileReadEndTime - fileReadStartTime);
                System.out.println(time + " " + fitnessGenBest + " " + fitnessRunBest + " " + runBESTTotalOverlap + " " + runBESTnContig + " " + duration);
                
                if(numOfGenWithNoUpgrade > nSSIter)
                    break;
            }
            if(fitnessRunBest > fitnessBest){
                BEST = runBEST;
                fitnessBest = fitnessRunBest;
            }
        }
         //////////////////////////////////////////////////////////////////////////////////////////////
        //MAIN LOOP ENDS HERE
        /////////////////////////////////////////////////////////////////////////////////////////////
        if(BEST == null){
           System.out.println("Holy");
        } 
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(strResultFile, true));

            long nContig = Utility.contigCount(BEST, overlapArray, threshold);
            List<String> contigs = Utility.contigs(BEST, overlapArray, threshold, tree);
            //System.err.println("contigs list size " + contigs.size());
            //bw.append(runId + " " + fitnessRunBest + " " + totalOverlap + " " + nContig + " " + duration);
            //bw.newLine();

            //loop to print fragment sequence
            for(int l = 0;l<contigs.size();l++) {
                //giving name to each sequence of .fna output file
                String seqName = ">Sequence_"+Integer.toString(l);
                bw.append(seqName);
                bw.newLine();

                String sequence = contigs.get(l);
                //printing each sequence. 70 bp per line
                for(int l1 = 0;l1<sequence.length();l1+=70) {
                    String s;
                    if(l1+70<sequence.length()) {
                        s = sequence.substring(l1, l1+70);
                    }
                    else{
                        s = sequence.substring(l1);
                    }
                    bw.append(s);
                    bw.newLine();
                }
            }

            bw.close();
        } catch (Exception e) {
            //System.err.println("Failed to created/open output file");
            System.err.println(e);

        }
    }

    public static GeneralizedSuffixTree ReadAndStoreFragments(String strFragmentFile)
    {
        GeneralizedSuffixTree tree = new GeneralizedSuffixTree();
        String strFragment = "";
        int nOrigStr;

        try
        {
            BufferedReader br = new BufferedReader(new FileReader(strFragmentFile));
            for(String line; (line = br.readLine()) != null; ) {
                if (line.charAt(0) == '>')
                {
                    if (!strFragment.isEmpty())
                    {
                        tree.put(strFragment);
                        strFragment = "";
                    }
                }
                else
                {
                    strFragment += line;
                }
            }
            tree.put(strFragment);

            br.close();
        }
        catch (Exception e)
        {
            System.out.println(e);
        }

        //
        // Now add the RC strings
        //
        nOrigStr = tree.getStringCount();
        for (int i = 0; i < nOrigStr; i++)
        {
            String strOrg = tree.getString(i);
            StringBuilder strRC = (new StringBuilder(strOrg)).reverse();
            for (int j = 0; j < strRC.length(); j++)
            {
                char rc = '\0';
                switch (strRC.charAt(j))
                {
                    case 'A':
                        rc = 'T';
                        break;
                    case 'T':
                        rc = 'A';
                        break;
                    case 'C':
                        rc = 'G';
                        break;
                    case 'G':
                        rc = 'C';
                        break;
                }

                strRC.setCharAt(j, rc);
            }

            tree.put(strRC.toString());
        }
        return tree;
    }

    public static int[][] generatePopulation(int popoulationSize, int individualSize) {

        int[][] population = new int[popoulationSize][];
        for (int i = 0; i < popoulationSize; i++) {
            population[i] = Utility.generateIndividual(individualSize);
        }
        return population;
    }

    public static int[][] generateSeededPopulation(int popoulationSize, int individualSize, int[][] ov, long maxIter, long threshold) {

        int[][] population = new int[popoulationSize][];
        for (int i = 0; i < popoulationSize; i++) {
            population[i] = LocalSearch.HillClimbing(ov, Utility.generateIndividual(individualSize), maxIter, threshold);
        }
        return population;
    }

    public static int[][] generateDiversePopulation(int diversePopulationSize, int individualSize, int[][] seed) {

        double maxDiversity = 0;
        int seedSize = seed.length;

        for (int i = 0; i < seed.length; i++) {
            double diversity = Utility.diversity(seed, seed[i], seedSize);

            if (diversity > maxDiversity) {
                maxDiversity = diversity;
            }
        }

        maxDiversity /= 2;

        int[][] population = new int[diversePopulationSize + seed.length][];

        for (int i = 0; i < seedSize; i++) {
            population[i] = seed[i];
        }

        for (int i = 0; i < diversePopulationSize; i++) {

            double diversity = 0;
            int[] diverseIndividual;

            do {
                diverseIndividual = Utility.generateIndividual(individualSize);
                diversity = Utility.diversity(population, diverseIndividual, seedSize + i);
            } while (diversity < maxDiversity);

            population[seedSize + i] = diverseIndividual;
        }

        return population;
    }

    private static int[][] getTournamentSelectedIndividuals(int [][] population, final double[] prop, int n)
    {
        int[][] selectedPopulation = new int[n][];
        for(int i = 0; i < n; i++)
        {
            selectedPopulation[i] = Utility.tournamentSelection(population, 7, prop);
        }

        return selectedPopulation;
    }

    private static int[][] getTopRankedIndividuals(int[][] population, final double[] prop, int n) {
        //
        // Inner class for index-sorting prop array in descending order
        // the indices array (declaration follows) holds the indices.
        //
        Comparator<Integer> comparator = new Comparator<Integer>() {
            public int compare(Integer arg0, Integer arg1) {

                Double a = new Double(prop[arg0]);
                Double b = new Double(prop[arg1]);

                int cmp = b.compareTo(a);
                if (cmp == 0)
                    return arg0.compareTo(arg1);

                return cmp;
            }
        };
        
        int[][] selectedPopulation = new int[n][];
        Integer[] indices = new Integer[prop.length];
        //PriorityQueue pqIndices = new PriorityQueue(prop.length, comparator);
        for (int i = 0; i < indices.length; i++) {
            indices[i] = i;
            //pqIndices.add(i);
        }
        Arrays.sort(indices, comparator);

        for (int i = 0; i < n; i++) {
            selectedPopulation[i] = population[indices[i]];
            //selectedPopulation[i] = population[pqIndices.poll()];
        }

        return selectedPopulation;
    }

    public static int[][] getFitIndividuals(int[][] population, double[] fitnessArray, int n)
    {

        return getTopRankedIndividuals(population, fitnessArray, n);
        //return getTournamentSelectedIndividuals(population, fitnessArray, n);
    }

    public static int[][] getDiverseIndividuals(int[][] population, final double[] diversityArray, int m)
    {
        return getTopRankedIndividuals(population, diversityArray, m);
    }
}
