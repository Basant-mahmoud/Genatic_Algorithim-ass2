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

// Function mean square error
double calculateError(const vector<double>& coefficients, const vector<pair<double, double>>& data)
{
    double error = 0.0;
    for (const auto& point : data)
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

// Function to perform two-point crossover
void twoPointCrossover(const vector<double>& parent1, const vector<double>& parent2, vector<double>& child1, vector<double>& child2)
{
    int crossoverPoint1 = rand() % parent1.size();
    int crossoverPoint2 = rand() % parent1.size();
    if (crossoverPoint1 > crossoverPoint2)
    {
        swap(crossoverPoint1, crossoverPoint2);
    }

    for (int i = 0; i < crossoverPoint1; ++i)
    {
        child1[i] = parent1[i];
        child2[i] = parent2[i];
    }
    for (int i = crossoverPoint1; i < crossoverPoint2; ++i)
    {
        child1[i] = parent2[i];
        child2[i] = parent1[i];
    }
    for (int i = crossoverPoint2; i < parent2.size(); ++i)
    {
        child1[i] = parent1[i];
        child2[i] = parent2[i];
    }
}

// Function nonuniform mutation
void nonUniformMutate(vector<double>& chromosome, int generation)
{
    for (double& gene : chromosome)
    {
        double mutationRate = static_cast<double>(rand()) / RAND_MAX;
        if (mutationRate < pm)
        {
            double deltaLower = gene - LB;
            double deltaUpper = UB - gene;

            double r1 = static_cast<double>(rand()) / RAND_MAX;
            double y = (r1 <= 0.5) ? deltaLower : deltaUpper;

            double r = static_cast<double>(rand()) / RAND_MAX;
            double mutation = y * (1 - pow(r, pow((1.0 - generation / max_generations), b)));
            if (y == deltaLower)
            {
                gene -= mutation;
            }
            else
            {
                gene += mutation;
            }
        }
    }
}

// Function tournament selection
vector<double> tournamentSelection(const vector<vector<double>>& population, const vector<double>& fitness)
{
    vector<double> bestChromosome;
    double bestFitness = 0.0;
    for (int i = 0; i < tournament_size; ++i)
    {
        int candidateIndex = rand() % population_size;
        double candidateFitness = fitness[candidateIndex];
        if (candidateFitness > bestFitness)
        {
            bestFitness = candidateFitness;
            bestChromosome = population[candidateIndex];
        }
    }

    return bestChromosome;
}

// Function elitist replacement
void elitistReplacement(vector<vector<double>>& oldPopulation, vector<vector<double>>& newPopulation,
                        const vector<double>& oldFitness)
{
    // Find the index of the best chromosome in the old population
    //new population have chromosome that apply on it crossover and mutation so we take some chromosome from the old population
    //to complete the function replacement like sildes k chromosome and num of chromosome population
    int bestIndex = distance(oldFitness.begin(), max_element(oldFitness.begin(), oldFitness.end()));

    // Copy the best chromosome to the next generation without crossover or mutation
    newPopulation.push_back(oldPopulation[bestIndex]);

    // Specify the number of additional best chromosomes to keep
    int numElites = static_cast<int>(population_size * 0.2) - 1; // Assuming 20% elitism ratio

    // Create a copy of oldFitness to avoid modifying the original vector
    vector<double> tempOldFitness = oldFitness;

    // Find indices of additional best chromosomes using a loop
    vector<int> eliteIndices(numElites);
    for (int i = 0; i < numElites; ++i)
    {
        auto maxIt = max_element(tempOldFitness.begin(), tempOldFitness.end());
        eliteIndices[i] = distance(tempOldFitness.begin(), maxIt);
        *maxIt = -1.0; // Mark as visited in the temporary copy
    }

    // Copy the additional best chromosomes to the next generation without crossover or mutation
    for (int i = 0; i < numElites; ++i)
    {
        newPopulation.push_back(oldPopulation[eliteIndices[i]]);
    }
}
// Function genetic algorithm
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
        for (int i = 0; i < population_size; ++i)
        {
            newPopulation[i] = tournamentSelection(population, fitness);
        }

        // Crossover
        for (int i = 0; i < population_size; i += 2)
        {
            if (static_cast<double>(rand()) / RAND_MAX < pc)
            {
                twoPointCrossover(newPopulation[i], newPopulation[i + 1], newPopulation[i], newPopulation[i + 1]);
            }
        }

        // Mutation
        for (vector<double>& chromosome : newPopulation)
        {
            nonUniformMutate(chromosome, generation);
        }
        // Elitist Replacement
        elitistReplacement(population, newPopulation, fitness);

        population = newPopulation;
    }

    // Find the best solution
    double bestFitness = 0.0;
    vector<double> bestChromosome;
    for (int i = 0; i < population_size; ++i)
    {
        double currentFitness = 1.0 / calculateError(population[i], data);
        if (currentFitness > bestFitness)
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
    if (!inputFile.is_open())
    {
        cerr << "Error: Could not open the input file.\n";
        return 1;
    }
    ofstream outputFile("output.txt");
    if (!outputFile.is_open())
    {
        cerr << "Error: Could not open the output file.\n";
        return 1;
    }
    int numDatasets;
    inputFile >> numDatasets;
    for (int datasetIndex = 1; datasetIndex <= numDatasets; ++datasetIndex)
    {
        int numPoints, degree;
        inputFile >> numPoints >> degree;

        vector<pair<double, double>> data(numPoints);
        for (auto& point : data)
        {
            inputFile >> point.first >> point.second;
        }
        vector<double> bestCoefficients = geneticAlgorithm(data, degree);

        // Output
        outputFile << "Dataset " << datasetIndex << ":\n";
        outputFile << "Coefficients: ";
        for (double coefficient : bestCoefficients)
        {
            outputFile << coefficient << ' ';
        }
        outputFile << "\nMean Square Error: " << calculateError(bestCoefficients, data) << "\n\n";
    }

    inputFile.close();
    outputFile.close();

    return 0;
}