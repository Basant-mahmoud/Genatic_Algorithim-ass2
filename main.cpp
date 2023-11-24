
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>
using namespace std;
const double LB = -10.0;
const double UB = 10.0;
const int population_size = 100;
const int max_generations = 100;
const double pc = 0.7;
const double pm = 0.05;
const double b = 1.5;
const int tournament_size = 5;
const int eliteSize=0.2 * population_size;
//function mean square error
double calculateError(const vector<double>& coefficients, const vector<pair<double, double>>& data)
{
    double error = 0.0;
    for(const auto& point : data)
    {
        double predicted = coefficients[0] + coefficients[1] * point.first + coefficients[2] * pow(point.first, 2);
        error += pow(predicted - point.second, 2);
    }
    return error / data.size();
}

// Function initialize population
vector<double> initializeChromosome(int degree)
{
    vector<double> chromosome(degree + 1);
    for (double& gene : chromosome)
    {
        gene = LB + static_cast<double>(rand()) / RAND_MAX * (UB - LB);
    }
    return chromosome;
}

//Function to perform two-point crossover
void twoPointCrossover(const vector<double>& parent1, const vector<double>& parent2, vector<double>& child1,vector<double>& child2)
{
    int crossoverPoint1 = rand() % parent1.size();
    int crossoverPoint2 = rand() % parent1.size();
    if(crossoverPoint1 > crossoverPoint2)
    {
        swap(crossoverPoint1, crossoverPoint2);
    }

    for(int i = 0; i < crossoverPoint1; ++i)
    {
        child1[i] = parent1[i];
    }
    for(int i = crossoverPoint1; i < crossoverPoint2; ++i)
    {
        child1[i] = parent2[i];
    }
    for(int i = crossoverPoint2; i < parent2.size(); ++i)
    {
        child1[i] = parent1[i];
    }
    ///////////////////
    for(int i = 0; i < crossoverPoint1; ++i)
    {
        child2[i] = parent2[i];
    }
    for(int i = crossoverPoint1; i < crossoverPoint2; ++i)
    {
        child2[i] = parent1[i];
    }
    for(int i = crossoverPoint2; i < parent1.size(); ++i)
    {
        child1[i] = parent2[i];
    }
}
//Function nonuniform mutation
void mutation(vector<double>& individual,int generationNumber)
{
    double deltalower,deltauber;
    for (int j = 0; j < individual.size(); j++)
    {
        double y,mutation;
        // random number for check pm
        if ((double)rand() / RAND_MAX < pm)
        {
            // random number from 0 to 1 for r1
            double r1=(double)rand() / RAND_MAX;

            if(r1<=0.5)
            {
                deltalower=individual[j]- LB;
                y= deltalower;
                double r2=(double)rand() / deltalower;
                mutation=y*(1-pow(r2,pow((1-generationNumber/max_generations),b)));
                individual[j] -=  mutation;
            }
            else
            {
                deltauber=UB-individual[j];
                y= deltauber;
                double r2=(double)rand() / deltauber;
                mutation=y*(1-pow(r2,pow((1-generationNumber/max_generations),b)));
                individual[j] +=  mutation;
            }
        }
    }
}


//Function tournament selection
vector<double> tournamentSelection(const vector<vector<double>>& population, const vector<double>& fitness)
{
    vector<double> bestChromosome;
    double bestFitness = 0.0;
    for(int i = 0; i < tournament_size; ++i)
    {
        int candidateIndex = rand() % population_size;
        double candidateFitness = fitness[candidateIndex];
        if(candidateFitness > bestFitness)
        {
            bestFitness = candidateFitness;
            bestChromosome = population[candidateIndex];
        }
    }

    return bestChromosome;
}

// Function elitist replacement
void elitistReplacement(const vector<vector<double>>& oldPopulation, vector<vector<double>>& newPopulation,
                        const vector<double>& oldFitness, const vector<double>& newFitness)
{
    //best fitness in old fitness
    int bestIndex = distance(oldFitness.begin(), max_element(oldFitness.begin(), oldFitness.end()));
    double bestFitness = oldFitness[bestIndex];

    if(*max_element(newFitness.begin(), newFitness.end()) > bestFitness)
    {
        //replace worest in fitness with best in fitness
        int worstIndex = distance(newFitness.begin(), min_element(newFitness.begin(), newFitness.end()));
        newPopulation[worstIndex] = oldPopulation[bestIndex];
    }

}

//Function genetic algorithm
vector<double> geneticAlgorithm(const vector<pair<double, double>>& data, int degree)
{
    srand(time(0));

    vector<vector<double>> population(population_size);
    for (vector<double>& chromosome : population)
    {
        chromosome = initializeChromosome(degree);
    }

    for (int generation = 0; generation < max_generations; ++generation)
    {
        // Evaluate fitness
        vector<double> fitness(population_size);
        for (int i = 0; i < population_size; ++i)
        {
            fitness[i] = 1.0 / calculateError(population[i], data);
        }

        // Tournament Selection
        vector<vector<double>> newPopulation(population_size);
        for(int i = 0; i < population_size; ++i)
        {
            newPopulation[i] = tournamentSelection(population, fitness);
        }

        // Crossover
        for(int i = 0; i < population_size; i += 2)
        {
            if(static_cast<double>(rand()) / RAND_MAX < pc)
            {
                twoPointCrossover(newPopulation[i], newPopulation[i + 1], newPopulation[i],newPopulation[i + 1]);
            }
        }

        // Mutation
        for(vector<double>& chromosome : newPopulation)
        {
            mutation(chromosome, generation);
        }
        // Elitist Replacement
        elitistReplacement(population, newPopulation, fitness, vector<double>(population_size, 0.0));

        population = newPopulation;
    }

    // Find the best solution
    double bestFitness = 0.0;
    vector<double> bestChromosome;
    for(int i = 0; i < population_size; ++i)
    {
        double currentFitness = 1.0 / calculateError(population[i], data);
        if(currentFitness > bestFitness)
        {
            bestFitness = currentFitness;
            bestChromosome = population[i];
        }
    }

    return bestChromosome;
}

int main()
{
    ifstream inputFile("curve_fitting_input.txt");
    if(!inputFile.is_open())
    {
        cerr << "Error: Could not open the input file.\n";
        return 1;
    }
    ofstream outputFile("output.txt");
    if(!outputFile.is_open())
    {
        cerr << "Error: Could not open the output file.\n";
        return 1;
    }
    int numDatasets;
    inputFile>>numDatasets;
    for(int datasetIndex = 1; datasetIndex <= numDatasets; ++datasetIndex)
    {
        int numPoints, degree;
        inputFile>>numPoints>>degree;

        vector<pair<double, double>> data(numPoints);
        for(auto& point : data)
        {
            inputFile>>point.first>>point.second;
        }
        vector<double> bestCoefficients = geneticAlgorithm(data, degree);

        // Output
        outputFile << "Dataset " << datasetIndex << ":\n";
        outputFile << "Coefficients: ";
        for(double coefficient : bestCoefficients)
        {
            outputFile << coefficient << ' ';
        }
        outputFile << "\nMean Square Error: " << calculateError(bestCoefficients, data) << "\n\n";
    }

    inputFile.close();
    outputFile.close();

    return 0;
}