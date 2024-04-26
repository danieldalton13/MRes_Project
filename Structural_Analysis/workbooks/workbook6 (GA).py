import numpy as np
import random
import math
import TrussSolver as ts  # This assumes your TrussSolver module is named TrussSolver.py
import pygame

class TrussGA:
    def __init__(self, nodes, loads, num_individuals=100, max_generations=100):
        self.nodes = nodes
        self.loads = loads
        self.num_individuals = num_individuals
        self.max_generations = max_generations
        self.population = [self.create_random_individual() for _ in range(num_individuals)]
        self.best_solution = None
        self.best_fitness = float('inf')

    def create_random_individual(self):
        # Create a random individual based on possible node connections
        individual = []
        for i in range(len(self.nodes)):
            # Each node can potentially connect to any other node, we restrict connections
            individual.append([random.choice([0, 1]) for _ in range(len(self.nodes))])
        return individual

    def fitness(self, individual):
        # Convert individual (genome) to connections compatible with the Truss class
        connections = []
        for i, connections_per_node in enumerate(individual):
            connected_nodes = [j for j, gene in enumerate(connections_per_node) if gene == 1 and i != j]
            connections.append(connected_nodes)

        truss = ts.Truss(self.nodes, connections)
        for load in self.loads:
            truss.addLoad(load[0], load[1])
        truss.solve()

        if truss.fmax > 18:
            return float('inf')  # Invalid truss configuration
        return truss.truss_weight()

    def select(self):
        # Tournament selection
        tournament_size = 5
        tournament = random.sample(self.population, tournament_size)
        tournament_fitness = [self.fitness(individual) for individual in tournament]
        return tournament[np.argmin(tournament_fitness)]

    def crossover(self, parent1, parent2):
        # Single-point crossover
        point = random.randint(1, len(parent1) - 1)
        child1 = parent1[:point] + parent2[point:]
        child2 = parent2[:point] + parent1[point:]
        return random.choice([child1, child2])

    def mutate(self, individual, mutation_rate=0.05):
        # Bit-flip mutation
        for i in range(len(individual)):
            for j in range(len(individual[i])):
                if random.random() < mutation_rate:
                    individual[i][j] = 1 - individual[i][j]
        return individual

    def run(self):
        for generation in range(self.max_generations):
            new_population = []
            while len(new_population) < self.num_individuals:
                parent1, parent2 = self.select(), self.select()
                offspring = self.crossover(parent1, parent2)
                offspring = self.mutate(offspring)
                new_population.append(offspring)
            self.population = new_population
            # Check best individual in current generation
            current_best = min(self.population, key=self.fitness)
            current_best_fitness = self.fitness(current_best)
            if current_best_fitness < self.best_fitness:
                self.best_fitness = current_best_fitness
                self.best_solution = current_best
            print(f"Generation {generation}: Best Fitness = {self.best_fitness}")

# Define nodes and loads
define_nodes = np.array([    
    [0,0],      #0
    [4,0],      #1
    [8,0],      #2
    [12,0],     #3
    [16,0],     #4
    [0,4],      #5
    [4,4],      #6
    [8,4],      #7
    [12,4],     #8
    [16,4],     #9
    [0,8],      #10
    [4,8],      #11
    [8,8],      #12
    [12,8],     #13
    [16,8],     #14
])

define_loads = [
    (0, [0, 25]),  
    (4, [0, 25]),  
    (12, [0, -50])
]

# Instantiate and run the GA
ga_truss = TrussGA(define_nodes, define_loads, num_individuals=100, max_generations=100)
ga_truss.run()
print("Optimal configuration found with weight:", ga_truss.best_fitness)

# At the end of your GA run, create a Truss object from the best solution
best_truss = ts.Truss(define_nodes, ga_truss.best_solution)
for load in define_loads:  # Applying each load to the truss
    best_truss.addLoad(load[0], load[1])
best_truss.solve()  # Solve the truss calculations if not already done

# Now visualize using PgTruss
truss_draw = ts.PgTruss(best_truss, 1600)
truss_draw.drawNodes()

# Keep the window open until the user closes it
running = True
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

pygame.quit()
