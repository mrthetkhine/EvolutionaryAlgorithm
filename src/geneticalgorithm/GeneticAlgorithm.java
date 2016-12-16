/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package geneticalgorithm;

import java.util.*;
/**
 *
 * @author student
 */
public class GeneticAlgorithm {


    int noOfParents ;
   
    int populationSize;
  

    ArrayList<Chromosome> populations = new ArrayList<Chromosome>(); //Store chromosome, chromosome are stored as integer array
    ArrayList<Chromosome> parents ;
    ArrayList<Chromosome> offsprings ;

    static double crossOverProab = 0.8;
    static double mutationProab = 0.3;
    
    Chromosome bestChromosome ;
    float bestFitness = -100000;
    GeneticAlgorithm(int populationSize)
    {
       
        this.populationSize = populationSize;
        this.noOfParents = 20;
    }

    void createChromosome()
    {
        for(int i=0; i< this.populationSize; i++)
        {
            int[] genes = new int[4];
            Chromosome chromosome = new Chromosome(genes);
            for(int g =0; g < genes.length;g++)
            {
              
                genes[g] = getRandom();
                
            }
            this.populations.add(chromosome);
          
        }
    }
    int getRandom()
    {
        Random r = new Random();
        return r.nextInt(30);
    }
    Chromosome getBestChromosome()
    {
        int best = 0;

        Chromosome bestCh = this.populations.get(0) ;
        boolean ok = false;
        
        for(int i=0;i< this.populations.size() ;i++)
        {
            Chromosome current = this.populations.get(i);
            if( current.fitness > bestCh.fitness)
            {
                bestCh = current;
            }
        }
  
        Chromosome ch = new Chromosome(bestCh.genes.length);
        
        for(int i=0;i< ch.genes.length;i++)
        {
            ch.genes[i]= bestCh.getGeneAt(i);
        }
        ch.evaluateFitness();
        return bestCh;

    }
    public void selectParent()
    {
        this.parents = new ArrayList<Chromosome>();

        double totalFitness = 0;
        for(int i=0;i < this.populations.size();i++)
        {
            totalFitness += this.populations.get(i).fitness;
        }
        double start = 0;

        for(int i=0;i < this.populations.size();i++)
        {
            Chromosome c = this.populations.get(i);
            double end = start + c.fitness / totalFitness;

            c.from = start;
            c.to = end;
            //System.out.println("Start "+c.from + " End "+c.to);
            start = end;

        }
        for(int i=0; i < this.noOfParents; i++)
        {
            double random = Math.random();
            //System.out.println("Selecting "+random);
            int p = 0;
            for(int c = 0; c< this.populations.size();c++)
            {
                Chromosome chrom = this.populations.get(c);
                if( chrom.from >= random &&  random < chrom.to)
                {
                    p = c;
                }
            }

            Chromosome best = this.populations.get(p);
            //System.out.println("Select Parent "+best);
            this.parents.add(best);

              
        }
    }
    public void evaluateFitnessForAllChromosome()
    {
        for(Chromosome c : this.populations)
        {
            c.evaluateFitness();
        }
        if(this.getBestChromosome().getFitness() > this.bestFitness)
        {
      
            this.bestChromosome = this.getBestChromosome();
            this.bestFitness = (float)this.getBestChromosome().getFitness();
        }
        
    }
    
    static Chromosome[] crossOver(Chromosome parentOne,Chromosome parentTwo)
    {
        Chromosome[] newChilds = new Chromosome[2];
        double crossOver = Math.random();

        int[] geneOne = new int[parentOne.getChromosomeLength()];
        int[] geneTwo = new int[parentTwo.getChromosomeLength()];
        //Copy gene into new chromosome
        for(int i=0;i< parentOne.getChromosomeLength();i++)
        {
            geneOne[ i ] = parentOne.getGeneAt(i);
            geneTwo[ i ] = parentTwo.getGeneAt(i);
        }
        Chromosome childOne = new Chromosome(geneOne);
        Chromosome childTwo = new Chromosome(geneTwo);

        if( crossOver < GeneticAlgorithm.crossOverProab )
        {
            int crossOverPoint = (int) Math.random() * (parentOne.getChromosomeLength() -1 );

            HashMap<Integer,Integer> mapping =new HashMap<Integer,Integer>();
            
            //Perform Corossover
            for(int c = crossOverPoint; c < parentOne.getChromosomeLength(); c++)
            {
                int gOne = parentOne.getGeneAt(c);
                int gTwo = parentTwo.getGeneAt(c);

                mapping.put(gOne, gTwo);
                mapping.put(gTwo,gOne);
                
                childOne.genes[ c ] = gTwo;
                childTwo.genes[ c ] = gOne;
            }
            //Pefrom validation
            ArrayList<Integer> childOneList = new ArrayList<Integer>();
            for(int i=0;i< childOne.getChromosomeLength();i++)
            {
                int gene = childOne.getGeneAt(i);
                if(childOneList.contains( gene ))
                {
                    childOne.genes[ i ] = mapping.get(gene);
                }
                else
                {
                    childOneList.add(gene);
                }
            }
            //Pefrom validation
            ArrayList<Integer> childTwoList = new ArrayList<Integer>();
            for(int i=0;i< childTwo.getChromosomeLength();i++)
            {
                int gene = childTwo.getGeneAt(i);
                if(childTwoList.contains( gene ))
                {
                    childTwo.genes[ i ] = mapping.get(gene);
                }
                else
                {
                    childTwoList.add(gene);
                }
            }
        }
        
        newChilds[ 0 ] = childOne;
        newChilds[ 1 ] = childTwo;
        
        return newChilds;
    }
    Chromosome mutation(Chromosome chro)
    {
        for(int i=0;i< chro.getChromosomeLength();i++)
        {
            double ran = Math.random();
            if(ran < GeneticAlgorithm.mutationProab)
            {
                ArrayList<Integer> allGene = new ArrayList<Integer>();
                for(int j=0;j< chro.getChromosomeLength();j++)
                {
                    int g = chro.getGeneAt(j);
                    allGene.add(g);
                }
                
                chro.genes[ i ] = getRandom();
            }
        }
        return chro;
    }
    
    Chromosome[] performCrossOverAndMutation()
    {
        Chromosome[] newChilds = new Chromosome[2];

        //System.out.println("Parent size is " + this.parents.size());

        int parentOneIndex = (int)(Math.random() * this.parents.size());
        int parentTwoIndex;
        //Pick up random index;
        do
        {
            //System.out.println("Picking up");
            parentTwoIndex = (int)(Math.random() * this.parents.size());
        }while( parentOneIndex == parentTwoIndex);

        Chromosome parentOne = this.parents.get(parentOneIndex);
        Chromosome parentTwo = this.parents.get(parentTwoIndex);
        newChilds = crossOver(parentOne, parentTwo);

        for(Chromosome c: newChilds)
        {
            this.mutation(c);
        }
        return newChilds;
    }
    double getTotalFitness()
    {
        double sum = 0;
        for(Chromosome c : this.populations)
        {
            sum += c.getFitness();
        }
        return sum;
    }
   
    public void runGA(int noOfIteration)
    {
        //this.createChromosome();

        for(int i=0;i< noOfIteration;i++)
        {
            this.evaluateFitnessForAllChromosome();
            this.selectParent();
            System.out.println("Iteration  " +i + " total fitness " + this.getTotalFitness() + " Best " +this.getBestChromosome().fitness);
            this.offsprings = new ArrayList<Chromosome>();
            for(int j=0;j < this.populationSize/2; j++ )
            {
                Chromosome [] childs = this.performCrossOverAndMutation();
                for(Chromosome c : childs)
                {
                    this.offsprings.add(c);
                }
            }
            this.populations = this.offsprings;
           
        }
    }
   
    void displayChromosome()
    {
        for(int i=0; i< this.populationSize; i++)
        {
            Chromosome chrom = this.populations.get(i);


            System.out.print("Chromosome " + i +" : "+ chrom);
            
            System.out.println();
            chrom.evaluateFitness();
            System.out.println("Fitness " + chrom.getFitness()  );
        }
    }
    
    public static void main(String[]args)
    {
        GeneticAlgorithm ga = new GeneticAlgorithm(40);
        ga.createChromosome();
        ga.runGA(200);
        Chromosome ch = ga.bestChromosome;
        
        int a = ch.genes[0], b = ch.genes[1], c = ch.genes[2], d = ch.genes[3];
        double result = a + 2*b + 3*c + 4*d - 30;
        System.out.println("Equation "+ result);
        System.out.println("A:"+a + " b:"+b + " c:"+ c + " d:" + d);
        
    }
}
