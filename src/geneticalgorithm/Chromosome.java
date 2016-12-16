/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package geneticalgorithm;

/**
 *
 * @author Mr T.Khine
 */
public class Chromosome {
    int[] genes;
   
    double fitness;
    
    double from,to; //To implement roulette wheel selection

    public Chromosome(int[] genes) 
    {
        this.genes = genes;
    }
    public Chromosome(int size)
    {
        this.genes = new int [size];
    }
    int getGeneAt(int index)
    {
        return this.genes[index];
    }
    double evaluateFitness()
    {
        //Equation to optimize is f  (x) =((a +  2b  +  3c  +  4d) -  30).
        double a = this.genes[0], b = this.genes[1], c = this.genes[2], d = this.genes[3];
        this.fitness = 1/( a + 2*b + 3*c + 4*d - 30);
        return this.fitness;
    }
    double getFitness()
    {
        return this.fitness;
    }
    int getChromosomeLength()
    {
        return this.genes.length;
    }
    @Override
    public String toString()
    {
        String s ="[";
        for(int i=0;i< this.genes.length;i++)
        {
            s+= this.genes[i] + " ";
        }
        s+= "]";
        return s;
    }
}
