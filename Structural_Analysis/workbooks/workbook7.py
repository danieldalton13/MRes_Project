import numpy as np
import random
import math
import TrussSolver as ts  # Assuming TrussSolver contains the necessary classes
import pygame

class TrussGA:
    def __init__(self, nodes, loads, connections, num_individuals=100, max_generations=100):
        self.nodes = nodes
        self.loads = loads
        self.connections = connections  # Add this line to initialize connections
        self.num_individuals = num_individuals
        self.max_generations = max_generations
        self.population = [self.create_random_individual() for _ in range(num_individuals)]
        self.best_solution = None
        self.best_fitness = float('inf')

    def create_random_individual(self):
        # Create a random individual based on specified connections
        individual = []
        for conn in self.connections:
            individual.append(random.choice([0, 1]))
        return individual

    def fitness(self, individual):
        # Convert individual to connections compatible with the Truss class
        connections = [[self.connections[i][j] for j, gene in enumerate(individual) if gene == 1] for i in range(len(self.connections))]
        truss = ts.Truss(self.nodes, connections)
        for load in self.loads:
            truss.addLoad(load[0], load[1])
        truss.solve()
        if truss.fmax > 18:
            return float('inf')  # Invalid truss configuration
        return truss.truss_weight()

    def select(self):
        tournament = random.sample(self.population, 5)
        return min(tournament, key=self.fitness)

    def crossover(self, parent1, parent2):
        point = random.randint(1, len(parent1) - 1)
        return random.choice([parent1[:point] + parent2[point:], parent2[:point] + parent1[point:]])

    def mutate(self, individual, mutation_rate=0.05):
        return [1 - gene if random.random() < mutation_rate else gene for gene in individual]

    def run(self):
        for generation in range(self.max_generations):
            new_population = []
            while len(new_population) < self.num_individuals:
                offspring = self.mutate(self.crossover(self.select(), self.select()))
                new_population.append(offspring)
            self.population = new_population
            current_best = min(self.population, key=self.fitness)
            current_best_fitness = self.fitness(current_best)
            if current_best_fitness < self.best_fitness:
                self.best_fitness = current_best_fitness
                self.best_solution = current_best
                print(f"New best fitness at Generation {generation}: {self.best_fitness}")
            else:
                print(f"Generation {generation}: Best Fitness = {self.best_fitness}")
        
        return self.best_solution

# Define nodes, loads, and connections
define_nodes = np.array([
    [0,0], [4,0], [8,0], [12,0], [16,0],
    [0,4], [4,4], [8,4], [12,4], [16,4],
    [0,8], [4,8], [8,8], [12,8], [16,8]
])

define_loads = [
    (0, [0, 25]), (4, [0, 25]), (12, [0, -50])
]

define_connections = [
    [1,5,6], [0,2,6], [1,3,6,7,8], [2,4,8], [3,8,9],
    [0,6,10], [0,1,2,5,7,10,11,12], [2,6,8,12], [2,3,4,7,9,12,12,14], [4,8,14],
    [5,6,11], [6,10,12], [6,7,8,11,13], [8,12,14], [8,9,13]
]

# Instantiate and run the GA
ga_truss = TrussGA(define_nodes, define_loads, define_connections, num_individuals=100, max_generations=100)
best_solution = ga_truss.run()

if best_solution:
    best_truss = ts.Truss(define_nodes, [[ga_truss.connections[i][j] for j, gene in enumerate(best_solution) if gene == 1] for i in range(len(best_solution))])
    for load in define_loads:
        best_truss.addLoad(load[0], load[1])
    best_truss.solve()
    
    pygame.init()
    truss_draw = ts.PgTruss(best_truss, 1600)
    truss_draw.drawNodes()
    
    running = True
    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
    
    pygame.quit()
else:
    print("No valid solution was found.")
