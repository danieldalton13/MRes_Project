import numpy as np
import random
import TrussSolver as ts  # Ensure TrussSolver includes necessary classes and methods

class TrussGA:
    def __init__(self, nodes, loads, connections, num_individuals=100, max_generations=100):
        self.nodes = nodes
        self.loads = loads
        self.connections = connections
        self.num_individuals = num_individuals
        self.max_generations = max_generations
        self.population = [self.create_random_individual() for _ in range(num_individuals)]
        self.best_solution = None
        self.best_fitness = float('inf')

    def create_random_individual(self):
        # Create a random individual by deciding for each connection if it is active (1) or not (0)
        individual = []
        for connection in self.connections:
            individual.append([random.choice([0, 1]) for _ in connection])
        return individual

    def fitness(self, individual):
        # Generate a list of active connections based on the individual's genes
        active_connections = []
        for i, node_connections in enumerate(self.connections):
            for j, conn in enumerate(node_connections):
                if individual[i][j] == 1:
                    active_connections.append([i, conn])
        
        try:
            truss = ts.Truss(self.nodes, active_connections)
            for load in self.loads:
                truss.addLoad(load[0], load[1])
            truss.solve()
        except Exception as e:
            print(f"Error during truss computation: {e}")
            return float('inf')  # Return a high penalty for unsolvable configurations

        if truss.fmax > 18:
            return float('inf')  # Invalid configuration with excessive force
        return truss.truss_weight()  # Minimize weight

    def select(self):
        tournament = random.sample(self.population, 5)
        return min(tournament, key=self.fitness)

    def crossover(self, parent1, parent2):
        # One-point crossover that respects node connection lists
        point = random.randint(0, len(self.connections) - 1)
        child = parent1[:point] + parent2[point:]
        return child

    def mutate(self, individual, mutation_rate=0.05):
        mutated = [[1-gene if random.random() < mutation_rate else gene for gene in node_connection] for node_connection in individual]
        return mutated

    def run(self):
        for generation in range(self.max_generations):
            new_population = []
            while len(new_population) < self.num_individuals:
                parent1, parent2 = self.select(), self.select()
                offspring = self.mutate(self.crossover(parent1, parent2))
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

# Define nodes, loads, and connections as you would have them set
nodes = np.array([
    [0,0], [4,0], [8,0], [12,0], [16,0],
    [0,4], [4,4], [8,4], [12,4], [16,4],
    [0,8], [4,8], [8,8], [12,8], [16,8]
])

loads = [
    (0, [0, 25]), (4, [0, 25]), (12, [0, -50])
]

connections = [
    [1,5,6], [0,2,6], [1,3,6,7,8], [2,4,8], [3,8,9],
    [0,6,10], [0,1,2,5,7,10,11,12], [2,6,8,12], [2,3,4,7,9,12,13,14], [4,8,14],
    [5,6,11], [6,10,12], [6,7,8,11,13], [8,12,14], [8,9,13]
]

ga_truss = TrussGA(nodes, loads, connections, num_individuals=50, max_generations=50)
best_solution = ga_truss.run()

if best_solution:
    best_truss = ts.Truss(nodes, [[i, conn] for i, node_connections in enumerate(best_solution) for conn in connections[i] if node_connections[conn] == 1])
    for load in loads:
        best_truss.addLoad(load[0], load[1])
    best_truss.solve()
    print("Optimized truss configuration found and solved.")
else:
    print("No valid solution was found.")
