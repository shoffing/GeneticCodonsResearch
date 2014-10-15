from numpy import *
import random
import math
import copy
import sys
import time
import Gnuplot, Gnuplot.funcutils

#=========
# Loading
#=========

# Information about the problem space that we're searching through
probInfo = {
    'order': [],
    'visits': {},
    'aminoForCodon': {},
    'codonsForAmino': {}
}

# Load the amino-codon pairings from the input file
with open('input_amino_codons', 'r') as f:
    for line in f:
        line = line.rstrip('\n')
        codon = line.split(',')[0]
        amino = line.split(',')[1]
        probInfo['aminoForCodon'][codon] = amino
        if amino in probInfo['codonsForAmino']:
            probInfo['codonsForAmino'][amino].append(codon)
        else:
            probInfo['codonsForAmino'][amino] = [codon]

# Load the codon link scores from the input file
scores = {}
with open('input_codon_scores', 'r') as f:
    for line in f:
        line = line.rstrip('\n')
        scores[line.split(',')[0]] = math.log(float(line.split(',')[2]) / float(line.split(',')[1]))

# Load the example DNA sequence, store initial sequence for comparison later
initSeq = []
with open('seq_polio', 'r') as f:
    next(f)
    for line in f:
        line = line.rstrip('\n')
        for cc in range(0, len(line), 3):
            codon = line[cc:cc+3]
            if codon in probInfo['aminoForCodon']:
                amino = probInfo['aminoForCodon'][codon]
                if amino != '*':
                    initSeq.append(codon)
                    probInfo['order'].append(amino)
                    probInfo['visits'][codon] = probInfo['visits'][codon] + 1 if (codon in probInfo['visits']) else 1
                else:
                    break


#===================
# Genetic Functions
#===================

# Generate a new random solution given a problem.
def generateRandomSolution(problem):
    # Make deep copy
    prob = copy.deepcopy(problem)

    solution = []
    ca = 0 # current amino order index
    while sum(prob['visits'].values()) > 0:
        curAmino = prob['order'][ca]

        # Find possible codons
        possibleCodons = []
        for pc in prob['codonsForAmino'][curAmino]:
            if pc in prob['visits'] and prob['visits'][pc] > 0:
                possibleCodons.append(pc)

        # Select the next codon to visit at this amino (ca) from the list of possible codons
        nextCodon = None
        if len(possibleCodons) > 0:
            nextCodon = possibleCodons[int(random.random() * len(possibleCodons))]
            prob['visits'][nextCodon] -= 1
        else:
            # All codons in this amino have no more remaining visits, just select a random one.
            nextCodon = prob['codonsForAmino'][curAmino][int(random.random() * len(prob['codonsForAmino'][curAmino]))];

        solution.append(nextCodon)

        if ca < len(prob['order']) - 1:
            ca += 1
        else:
            ca = 0

    return solution


# Randomly mutate part (amt) of a given solution
def mutate(solution, amt, problem):
    # Make deep copy
    sol = copy.deepcopy(solution)

    for i in range(int(amt * len(sol))):
        mutIndex = int(random.random() * len(sol))

        # Find possible nodes to swap with
        pSwaps = []
        for j in range(len(sol)):
            if j != mutIndex and problem['aminoForCodon'][sol[j]] == problem['aminoForCodon'][sol[mutIndex]]:
                pSwaps.append(j)

        if len(pSwaps) > 0:
            # Swap the current node with a random node selected from pSwaps
            swapIndex = pSwaps[int(random.random() * len(pSwaps))]
            sol[mutIndex], sol[swapIndex] = sol[swapIndex], sol[mutIndex]

    return sol


# Given two parent solutions a and b, the number of crossover points
# to use, and the problem space: generates a child.
def crossover(a, b, numPoints, problem):
    # Get a list of crosspoints to use
    crossPoints = []
    for i in range(1, numPoints + 1):
        crossPoints.append( round(i * len(a) / float(numPoints + 1)) )

    # Perform crossover
    child = []
    isA = random.random() < 0.5 # Keep track of whether we're taking pairs from A or B
    ccp = 0 # current crossover point
    for i in range(len(a)):
        child.append(a[i] if isA else b[i])

        if len(crossPoints) > ccp and i > crossPoints[ccp]:
            isA = not isA
            ccp += 1

    if problem != None:
        fix(child, problem)

    return child


# Given a solution that doesn't match the problem criteria,
# fix that solution so that it does match the criteria
def fix(solution, problem):
    # Make deep copy of problem
    prob = copy.deepcopy(problem)

    for i in range(len(solution)):
        codon = solution[i]
        amino = prob['aminoForCodon'][codon]

        if prob['visits'][codon] > 0:
            prob['visits'][codon] -= 1
        else:
            # Find an alternative for this one
            possibleCodons = []
            for pc in prob['codonsForAmino'][amino]:
                if pc in prob['visits'] and prob['visits'][pc] > 0:
                    possibleCodons.append(pc)

            newCodon = None
            if len(possibleCodons) > 0:
                newCodon = possibleCodons[int(random.random() * len(possibleCodons))]
                prob['visits'][codon] -= 1
            else:
                # Just pick a random codon then
                newCodon = prob['codonsForAmino'][amino][int(random.random() * len(prob['codonsForAmino'][amino]))];

            solution[i] = newCodon



# Calculate the score of a solution using the score matrix.
def calcScore(solution):
    score = 0
    for i in range(len(solution) - 1):
        pair = solution[i] + solution[i+1]

        if pair in scores:
            score += scores[pair]
        else:
            print 'PAIR NOT FOUND IN SCORE MATRIX: ' + pair

    return score


#==========
# Testing
#==========

NUM_GENS = 100000 # Generations
NUM_SOLS_PER_GEN = 10 # solutions per generation
NUM_CROSSPOINTS = 2 # crossover points
NUM_SELECTIONS = 3 # Top solutions survive each generation
NUM_CHILDREN = 5 # children are created from those
NUM_RANDOM = NUM_SOLS_PER_GEN - (NUM_SELECTIONS + NUM_CHILDREN) # random solutions fill in the gaps
MUTATION_AMOUNT = 0.01 # Mutate everything that we bring over a bit
NUM_PRESERVE = 1 # Prevent the algorithm from mutating the best-performing solution(s). While unrealistic, it prevents us from going backwards.


# Create a new generation given a previous one
def newGeneration(pool):
    sortedPool = sorted(pool, key=calcScore)
    
    # Selection - Bring in the best performers
    newGen = sortedPool[:NUM_SELECTIONS]

    # Crossover - Breed the top performers
    for i in range(NUM_CHILDREN):
        # Select 2 random parents
        ai = int(random.random() * len(newGen))
        bi = int(random.random() * len(newGen))
        if ai == bi:
            bi = (bi + 1) % len(newGen)

        a = newGen[ai]
        b = newGen[bi]

        newGen.append(crossover(a, b, NUM_CROSSPOINTS, probInfo))

    # Mutation
    newGen[NUM_PRESERVE:] = map(mutate, newGen[NUM_PRESERVE:], [MUTATION_AMOUNT] * len(newGen[NUM_PRESERVE:]), [probInfo] * len(newGen[NUM_PRESERVE:]))

    # Random new solutions to fill the gap
    newGen.extend([generateRandomSolution(probInfo) for x in range(NUM_RANDOM)])

    return newGen


#for tv in arange(0.00, 0.20, 0.01):
#    MUTATION_AMOUNT = tv

start = time.time() # Measure elapsed time

# Store the best solution, this will be the output.
bestSolution = None

# store average scores and best scores of each generation
averageScores = []
bestScores = []

# Generate initial population
pool = [generateRandomSolution(probInfo) for x in range(NUM_SOLS_PER_GEN)]

poolScores = map(calcScore, pool)
averageScores.append(sum(poolScores) / float(len(poolScores)))
bestScores.append(sorted(poolScores)[0])

for i in range(NUM_GENS):
    sys.stdout.write('\rComplete: %g%%' % round(100 * (i+1) / float(NUM_GENS), 2) + (' '*6))
    sys.stdout.flush()

    pool = newGeneration(pool)

    # Save best solution
    curBestSolution = sorted(pool, key=calcScore)[0]
    if bestSolution == None or calcScore(curBestSolution) < calcScore(bestSolution):
        bestSolution = curBestSolution

    poolScores = map(calcScore, pool)
    averageScores.append(sum(poolScores) / float(len(poolScores)))
    bestScores.append(sorted(poolScores)[0])

end = time.time()


#=========
# Output
#=========

print '\n' + str(NUM_GENS) + ' generations in ' + str(round(end - start, 2)) + ' seconds [' + str(round(NUM_GENS / (end-start), 1)) + ' generations per second]'

# Outputting to text file
with open('output/ga_[' + ','.join(map(str,[NUM_GENS,NUM_SOLS_PER_GEN,NUM_CROSSPOINTS,NUM_SELECTIONS,NUM_CHILDREN,NUM_RANDOM,MUTATION_AMOUNT,NUM_PRESERVE])) + '].txt', 'w') as f:
    f.write('Best sequence [' + str(calcScore(bestSolution)) + ']:\n' + ''.join(bestSolution) + '\n\n')
    f.write('Initial sequence [' + str(calcScore(initSeq)) + ']:\n' + ''.join(initSeq) + '\n\n')

    #f.write('Gen\tAvg\tBest\n')
    #for i in range(NUM_GENS + 1):
    #    f.write(str(i) + '\t' + str(averageScores[i]) + '\t' + str(bestScores[i]) + '\n')


# Graphing with Gnuplot
g = Gnuplot.Gnuplot(debug=0, persist=0)
g.title('GA Performance [' + ', '.join(map(str,[NUM_GENS,NUM_SOLS_PER_GEN,NUM_CROSSPOINTS,NUM_SELECTIONS,NUM_CHILDREN,NUM_RANDOM,MUTATION_AMOUNT,NUM_PRESERVE])) + ']')
g('set xlabel "Generation #"')
g('set ylabel "Score"')
g('set xrange [:' + str(NUM_GENS) + ']')
#g('set yrange [-400:-20]')
g('set logscale x')

g('set terminal png')

g('set style data steps')

best = Gnuplot.Data(zip(range(0, len(bestScores)), bestScores), title='Best Score')
avg = Gnuplot.Data(zip(range(0, len(averageScores)), averageScores), title='Average Score')

g.plot(best, avg)
g.hardcopy('output/ga_[' + ','.join(map(str,[NUM_GENS,NUM_SOLS_PER_GEN,NUM_CROSSPOINTS,NUM_SELECTIONS,NUM_CHILDREN,NUM_RANDOM,MUTATION_AMOUNT,NUM_PRESERVE])) + '].png', terminal='png')

g.close()


raw_input('Please press return to continue...\n')