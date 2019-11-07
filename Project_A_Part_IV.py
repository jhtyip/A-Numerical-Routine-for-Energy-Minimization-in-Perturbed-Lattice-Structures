# The Chinese University of Hong Kong, Department of Physics
# PHYS4061 Project A Part IV
# Written by Yip Hoi Tung (1155092406) in Python 2.7
# December 2018

import sys
sys.path.append("Project_A_Part_I") # Change path to the folder where the Part I python script is under
# Import modified (added the BC8 structure and other slight modifications) functions from Part I
import Project_A_Part_I_modified_LJ as pAp1_LJ # Tailor-made for LJ potential
import Project_A_Part_I_modified_Ter as pAp1_Ter # Tailor-made for Tersoff potential

from random import uniform # For the random perturbation
import math

from copy import deepcopy # For copying nested lists

#################### Functions from Part II (slightly modified) ####################
def dot(u, v): # Dot product of 2 vectors; each vector is a list
    dotSum = 0
    for i in range(len(u)):
        dotSum += u[i]*v[i]
    
    return dotSum

def cross(u, v): # Cross product of 2 vectors; each vector is a list
    return [u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0]]

def recipVect(a): # Calculate the reciprocal lattice vectors b from the lattice vectors a
    b = []
    b.append([i/dot(a[0], cross(a[1], a[2])) for i in cross(a[1], a[2])])
    b.append([i/dot(a[0], cross(a[1], a[2])) for i in cross(a[2], a[0])])
    b.append([i/dot(a[0], cross(a[1], a[2])) for i in cross(a[0], a[1])])

    return b # b is a list of vectors; each vector is a list

def coorTran(a, r): # Coordinate transform: transform the coordinates of a point r into fractional coordinates given lattice vectors a
    b = recipVect(a) # Get reciprocal lattice vectors b from lattice vectors a
    
    r_frac = [] # r_frac is the point r in fractional coordinates, i.e. [n1, n2, n3] in Step 2
    for i in range(3): # Loop the 3 dimensions
        frac = dot(r, b[i]) % 1 # (r dot bi) mod 1
        if frac >= 0.5: # If the resulting fractional coordinate is >= 0.5, translate it to (frac - 1); otherwise keep it
            r_frac.append(frac - 1)
        else:
            r_frac.append(frac)

    r_cart = [] # r_cart is the Cartesian coordinates corresponding to the fractional coordinates of the point, i.e. [x2, y2, z2] in Step 2
    for i in range(3): # Loop the 3 dimensions
        cart = 0
        for j in range(3): # Loop the 3 lattice vectors: x2 = x-component of (n1*a1 + n2*a2 + n3*a3)
            cart += r_frac[j]*a[j][i]
        r_cart.append(cart)

    return r_frac, r_cart # Return both the fractional and corresponding Cartesian coordinates of the point r

def infoFrmXyzFile(fileName): # This program takes data from a .xyz file which can be from Project A Part I. This function is to take useful data from the file.
    content = open(fileName, "r").readlines() # fileName is the name of the file, e.g. "inputFiles\Sc_3_3_3.xyz"
    
    comment = content[1] # The comment line

    atoms = [] # List of the atoms
    for i in range(2, len(content)): # Loop the whole atom list in the .xyz file
        if content[i][-1] == "\n": # This condition is due to a particular detail of the formatting, which is unimportant
            atoms.append(content[i][:-1].split("\t")) # As of this point, each atom is stored as a list [element, x, y, z]
        else:
            atoms.append(content[i].split("\t"))

    for i in range(len(atoms)): # Loop the atom list 
        for j in range(1, 4): # Change the coordinates from string to float type
            atoms[i][j] = float(atoms[i][j])
        atoms[i].append(i + 1) # Number each atom (as the atom's ID). So, as of this point, each atom is stored as a list [element, x, y, z, ID]

    return comment, atoms

def latticeVec(latticeInfo): # Return corresponding lattice vectors according to the file name. latticeInfo is a list in this format: [lattice constant, periodicity in x, y, z]
    return [[latticeInfo[0]*latticeInfo[1], 0.0, 0.0], [0.0, latticeInfo[0]*latticeInfo[2], 0.0], [0.0, 0.0, latticeInfo[0]*latticeInfo[3]]]

def getDist(a, refAtom, candAtom): # Calculate separation distance of the "candidate atom" from the "reference atom". a is the lattice vectors, refAtom is the reference atom and candAtom is the candidate atom.
    from math import sqrt

    relDist = [] # Relative distances = [(delta x), (delta y), (delta z)]
    for i in range(1, 4): # Recall each atom has the format [element, x, y, z, ID], so index 1, 2, 3 are the coordinates
        relDist.append(candAtom[i] - refAtom[i])
    
    candAtom_frac, candAtom_cart = coorTran(a, relDist) # Get the coordinate transformation of the relative distances
    
    sqDist = 0
    for i in range(3): # distance^2 = x^2 + y^2 + z^2
        sqDist += candAtom_cart[i]**2

    return sqrt(sqDist) # Return the distance

def takeDist(atom): # This is for getting the distance element in the atom list for sorting in the next function
    return atom[5]

def getNeighborList(latticeInfo, atoms, distCutoff): # Generate the list of neighbors. latticeInfo and atoms are from the .xyz file from Part I
    a = latticeVec(latticeInfo) # Get lattice vectors

    neighborList = []
    for i in range(len(atoms)): # Loop all atoms in the .xyz file
        neighborList_i = [] # The neighbor list (neighborList) is composed of lists of neighbors for each atom (neighborList_i)
        
        neighborList_i.append(atoms[i]) # The first item in neighborList_i is the atom which its neighbors are going to be put in this list. That means neighborList_i = [[atom_i], [1st neighbor of atom_i], [2nd neighbor of atom_i], ...]
        
        neighborList_i_dummy = [] # This dummy variable is almost an unsorted version of neighborList_i
        for j in range(len(atoms)): # Loop all atoms
            dist = getDist(a, atoms[i], atoms[j]) # Calculate separation distance of the atom from the atom in interest
            if dist <= distCutoff and dist != 0: # If the separation distance is smaller than or equal to the distance cutoff and is non-zero (zero implies the same atoms were compared), the atom is a neighbor
                dummyAtom = [k for k in atoms[j]] # Copy the atom. Recall that each atom is a list: [element, x, y, z, ID]
                dummyAtom.append(dist) # Add the separation distance to the list. Now each atom is a list: [element, x, y, z, ID, distance]
                neighborList_i_dummy.append(dummyAtom) # Add the atom to the unsorted dummy list
        
        neighborList_i_dummy_sorted = sorted(neighborList_i_dummy, key=takeDist) # Sort the dummy list according to the distance elements of the atoms
        for k in range(len(neighborList_i_dummy_sorted)): # Add the atoms in the sorted list to neighborList_i
            neighborList_i.append(neighborList_i_dummy_sorted[k])
        
        neighborList.append(neighborList_i) # Add neighborList_i to neighborList

    return neighborList

def maxNumOfNeighbors(neighborList): # Find the maximum number of neighbors for an atom in the neighbor list can have. All atoms should have the same number of neighbors in our highly symmetric case, so this is just for the generality
    numOfNeighbors = []
    for i in range(len(neighborList)):
        numOfNeighbors.append(len(neighborList[i]) - 1)

    return max(numOfNeighbors)

def getNeighborListFile(fileName, distCutoff, outputFileName): # Output to file. This function looks complicated but all are just the formatting.
    comment, atoms = infoFrmXyzFile(fileName) # Get info of the unit cell from .xyz file from Part I
    comment = comment.split("; ")

    outputFile = open(outputFileName, "w")

    outputFile.write(comment[0] + "\n") # Print info of the unit cell
    outputFile.write(comment[1] + "\n")
    outputFile.write(comment[2])
    outputFile.write("Total number of atoms: " + str(len(atoms)) + "\n\n" + "List of atoms:\n")
    
    outputFile.write("{: <10}{: <10}{: <10}{: <10}{: <10}".format("ID", "Element", "x-coor", "y-coor", "z-coor") + "\n") # Print list of atoms. The "{: <10}" strings control the width of the columns in the table
    for i in range(len(atoms)):
        outputFile.write("{: <10}{: <10}{: <10}{: <10}{: <10}".format(str(atoms[i][4]), str(atoms[i][0]), str(atoms[i][1]), str(atoms[i][2]), str(atoms[i][3])) + "\n")
    
    outputFile.write("\nNeighbor list (distance cutoff = %s Angstrom):\nAtom ID\t\tIDs of the atom's neighbors (separation distance in parentheses, from nearest to farthest)\n \t\t" % distCutoff) # Print neighbor list
    
    latticeInfo = [float(comment[1].split("Lattice constant: ")[1].split(" Angstrom")[0])]
    latticeInfo.extend([float(i) for i in comment[2].split("Periodicity (x * y * z): ")[1].split("\n")[0].split(" * ")]) # latticeInfo is now a list in this format: [lattice constant, periodicity in x, y, z]

    neighborList = getNeighborList(latticeInfo, atoms, distCutoff) # Get neighbor list
    
    rowForPrint = [] # Contains the elements to be printed in a row
    columnFormat = "" # Contains "{: <15}" strings for formatting the width of the columns in the table
    for n in range(maxNumOfNeighbors(neighborList)): # This is to print the the first row which is the labeling
        rowForPrint.append("N%g" % (n + 1))
        columnFormat += "{: <15}"
    outputFile.write(columnFormat.format(*rowForPrint) + "\n")
    
    for i in range(len(neighborList)): # This is to print the entire neighbor list
        outputFile.write(str(neighborList[i][0][4]) + "\t\t")
        rowForPrint = []
        for j in range(1, len(neighborList[i])):
            rowForPrint.append(str(neighborList[i][j][4]) + "(" + str(round(neighborList[i][j][5], 2)) + ")")
        columnFormat = ""
        for j in range(len(rowForPrint)):
            columnFormat += "{: <15}"
        outputFile.write(columnFormat.format(*rowForPrint) + "\n")
        
    outputFile.close()
    print "The neighbor list (%s) has been generated." % outputFileName

    return neighborList, latticeInfo

############################# Functions from Part III ##############################
########################### For Lennard-Jones potential ############################
def get_U_LJ(neighborList): # LJ stands for Lennard-Jones, obviously
    #print "\nImplementing the Lennard-Jones potential:"
    
    # Parameters for Xenon
    epsilon = 320e-16
    sigma = 3.98
    
    # Parameters from inputs
    #epsilon = float(raw_input("Input parameter epsilon (in erg): ")) # Parameters
    #sigma = float(raw_input("Input parameter sigma (in Angstrom): "))

    U_tot_LJ = 0.0 # Total potential energy of the lattice
    for i in range(len(neighborList)): # Uncomment this line and the line below when the lattice structure is not symmetric enough, i.e. energies for each atom are not all the same
    #for i in range(1): # Equivalent to setting i = 0, i.e. Look at the neighborList[0] which is the neighbor list of the first atom
        distList = [neighborList[i][x][5] for x in range(1, len(neighborList[i]))] # Extract the distances of all the neighbors and make a list. e.g. neighborList[0][2][5]: In the neighbor list of atom 0, the 2nd neighbor (index 0 is the atom 0 itself) which has the separation distance stored in index 5
        for j in range(len(distList)): # Loop the distances
            U_tot_LJ += (sigma/distList[j])**12 - (sigma/distList[j])**6 # Summation of the LJ potential
    U_tot_LJ *= 0.5*4.0*epsilon
    #U_tot_LJ *= 0.5*len(neighborList)*4.0*epsilon # 0.5 for compensating the double counting; len(neighborList) for the number of atoms in the lattice
    U_avg_LJ = U_tot_LJ/len(neighborList) # Energy per atom

    # Uncomment the below to calculate constants A and B (see report attached)
    '''
    minDist = min(distList) # Nearest neighbor's distance
    distList = [i/minDist for i in distList]
    checkSumA = 0.0 # A
    checkSumB = 0.0 # B
    for i in range(len(distList)):
        checkSumA += distList[i]**-12
        checkSumB += distList[i]**-6
    
    print "A is %g" % checkSumA
    print "B is %g" % checkSumB
    '''

    return U_tot_LJ, U_avg_LJ # Return total energy and energy per atom

############################## For Tersoff potential ###############################
def getLatticeInfo(fileName): # This function reads the .xyz file from Part I and get information of the lattice and output a list. The ultimate reason for this is to compose the lattice vector, which is done by calling the function latticeVec(latticeInfo) which takes this function's output as input
    comment, atoms = infoFrmXyzFile(fileName) # Call a function written in Part II to read the .xyz file
    comment = comment.split("; ")

    latticeInfo = [float(comment[1].split("Lattice constant: ")[1].split(" Angstrom")[0])]
    latticeInfo.extend([float(i) for i in comment[2].split("Periodicity (x * y * z): ")[1].split("\n")[0].split(" * ")]) # latticeInfo is now a list in this format: [lattice constant, periodicity in x, y, z]

    return latticeInfo

def get_U_Ter(neighborList, fileName): # Ter stands for Tersoff, obviously
    U_tot_Ter = 0.0 # Total energy of the lattice
    for i in range(len(neighborList)): # Uncomment this line and the line below when the lattice structure is not symmetric enough, i.e. energies for each atom are not all the same
    #for i in range(1): # Equivalent to setting i = 0, i.e. Look at the neighborList[0] which is the neighbor list of the first atom
        neighborList_i = neighborList[i] # neighborList_i is a list of neighbor atoms of atom i
        for j in range(1, len(neighborList_i)): # Loop the neighbor atoms
            U_tot_Ter += V_ij(neighborList_i, j, fileName) # Add potential energy of each pair (i & j) to the total energy
    U_tot_Ter *= 0.5
    #U_tot_Ter *= 0.5*len(neighborList) # 0.5 for compensating the double counting; len(neighborList) for the number of atoms in the lattice
    U_avg_Ter = U_tot_Ter/len(neighborList) # Energy per atom

    return U_tot_Ter, U_avg_Ter # Return total energy and energy per atom

# The functions below are named the same as those used by J. Tersoff in his paper (see report attached) 

def V_ij(neighborList_i, j, fileName): # Tersoff potential energy due to the pair atom i & j; j is one of the neighbors of i
    r = neighborList_i[j][5] # Atom j, index 5 (which stores the separation distance)

    if f_c(r) == 0.0: # If the cutoff function returns zero (i.e. cut off), there is no need to calculate other terms
        return 0.0
    else:
        return f_c(r)*(a_ij(neighborList_i, j)*f_R(r) + b_ij(neighborList_i, j, fileName)*f_A(r)) # The Tersoff potential

def f_c(r): # Cutoff function
    from math import sin, pi

    if r <= R - D:
        return 1.0
    elif r < R + D and r > R - D: # So, the smooth cutoff takes place from (R - D) to (R + D)
        return 0.5 - 0.5*sin(pi/2.0*(r - R)/D)
    elif r >= R + D:
        return 0.0

def f_R(r):
    from math import exp

    return A*exp(-1.0*lambda_1*r)

def f_A(r):
    from math import exp

    return -1.0*B*exp(-1.0*lambda_2*r)

def b_ij(neighborList_i, j, fileName):
    from math import exp

    zeta_ij = 0.0
    for k in range(1, len(neighborList_i)): # Loop the neighbor list for the third atom (k) is introduced
        if k != j: # Atom k cannot be the same as j 
            zeta_ij += f_c(neighborList_i[k][5])*g(neighborList_i[0], neighborList_i[j], neighborList_i[k], fileName)*exp((lambda_3**3)*((neighborList_i[j][5] - neighborList_i[k][5])**3))

    return (1.0 + (beta**n)*(zeta_ij**n))**(-1.0/2.0/n)

def g(atom_i, atom_j, atom_k, fileName): # g is a function of theta in the literature. To calculate theta here, the coordinates of the atom i, j & k are used
    vector_ij = [atom_j[x] - atom_i[x] for x in range(1, 4)] # Relative position
    vector_ik = [atom_k[x] - atom_i[x] for x in range(1, 4)]

    vector_ij_frac, vector_ij_cart = coorTran(latticeVec(getLatticeInfo(fileName)), vector_ij) # Center at atom i, get the coordinates of j & k after imposing the periodic boundary condition
    vector_ik_frac, vector_ik_cart = coorTran(latticeVec(getLatticeInfo(fileName)), vector_ik)

    cosAngle = dot(vector_ij_cart, vector_ik_cart)/vectorMag(vector_ij_cart)/vectorMag(vector_ik_cart) # cos(theta) = (dot product)/[(length1)(length2)]

    return 1.0 + (c/d)**2 - c**2/(d**2 + (h - cosAngle)**2)

def vectorMag(vector): # To calculate the length of a vector
    from math import sqrt

    sumOfSquares = 0
    for i in range(len(vector)):
        sumOfSquares += vector[i]**2

    return sqrt(sumOfSquares)

def a_ij(neighborList_i, j):
    from math import exp

    eta_ij = 0.0
    for k in range(1, len(neighborList_i)): # Similar to and simplier than b_ij
        if k != j:
            eta_ij += f_c(neighborList_i[k][5])*exp((lambda_3**3)*((neighborList_i[j][5] - neighborList_i[k][5])**3))

    return (1.0 + (alpha**n)*(eta_ij**n))**(-1.0/2.0/n)

# Parameters for the Tersoff potential (unit in comment) (J. Tersoff published a few sets of parameters, this is the one updated in May 1988)
A = 1.8308e3 # eV
B = 4.7118e2 # eV
lambda_1 = 2.4799 # A**-1
lambda_2 = 1.7322 # A**-1
alpha = 0.0
beta = 1.0999e-6
n = 7.8734e-1
c = 1.0039e5
d = 1.6218e1
h = -5.9826e-1
lambda_3 = 1.7322 # A**-1
R = 2.85 # A
D = 0.15 # A

# An older set of parameters
# A = 3.2647e3
# B = 9.5373e1
# lambda_1 = 3.2394
# lambda_2 = 1.3258
# alpha = 0.0
# beta = 3.3675e-1
# n = 2.2956e1
# c = 4.8381
# d = 2.0417
# h = 0.0
# lambda_3 = 1.3258
# R = 3.0
# D = 0.2

############################## Functions for Part IV ###############################
def formatNeighborList(neighborList, latticeInfo): # Change coordinates of neighbors to coordinates by periodic boundary condition
    for i in range(len(neighborList)): # Loop all atoms' neighbor list
        for j in range(1, len(neighborList[i])): # Loop all neighbors in that list
            disVector = [neighborList[i][j][k] - neighborList[i][0][k] for k in range(1, 4)] # Displacement vector to that neighbor atom
            disVector_frac, disVector_cart = coorTran(latticeVec(latticeInfo), disVector) # Fractional coordinates and the corresponding cartisian coordinates of that displacement vector

            for k in range(1, 4):
                neighborList[i][j][k] = neighborList[i][0][k] + disVector_cart[k - 1] # Correct coordinates of each neighbor atom to the atom which it is the neighbor of

    return neighborList

def updateAtomCoor(neighborList, atomID, coorChange): # Update all the coordinates of an atom in the neighborList. The atom can have more than one coordinates because it may be a neighbor of another atom after imposing the periodic boundary condition
    for i in range(len(neighborList)):
        for j in range(len(neighborList[i])):
            if neighborList[i][j][4] == atomID: # Loop all atoms in neighborList, if it is the atom we are looking for
                for k in range(1, 4):
                    neighborList[i][j][k] += coorChange[k - 1] # Add the change to each coordinate

    return neighborList

def updateAtomDist(neighborList): # Update to the correct separation distances to neighbor atoms since their coordinates had changed
    for i in range(len(neighborList)):
        for j in range(1, len(neighborList[i])):
            neighborList[i][j][5] = math.sqrt(sum([(neighborList[i][j][k] - neighborList[i][0][k])**2 for k in range(1, 4)])) # Calculate the correct distance to neighbor atom since the coordinates of the neighbor atom had been updated
    return neighborList

def getPerturbedList(neighborList, latticeInfo, maxPerturbDist): # Generate the neighbor list after random perturbation of the lattice
    neighborList_stable = deepcopy(neighborList) # Make a copy of the original neighborList
    
    neighborList_perturbed = deepcopy(neighborList)
    for i in range(len(neighborList)): # Loop all atoms in the lattice
        atomID = neighborList[i][0][4] # Atom being perturbed at this instance
        perturb_x = uniform(-maxPerturbDist/100.0*latticeInfo[0], maxPerturbDist/100.0*latticeInfo[0]) # maxPerturbDist is in the unit "percentage of lattice constant". The uniform() function then randomly chooses (in a uniform distribution) a value between the negative and the positive maximum distance as the perturbation in a direction to be added to the atom's x, y, or z-coordinate
        perturb_y = uniform(-maxPerturbDist/100.0*latticeInfo[0], maxPerturbDist/100.0*latticeInfo[0])
        perturb_z = uniform(-maxPerturbDist/100.0*latticeInfo[0], maxPerturbDist/100.0*latticeInfo[0])

        coorChange = [perturb_x, perturb_y, perturb_z]
        neighborList_perturbed = updateAtomCoor(neighborList_perturbed, atomID, coorChange) # Update all the coordinates associated with this atom

    neighborList_perturbed = updateAtomDist(neighborList_perturbed) # Update the separation distances

    return neighborList_perturbed, neighborList_stable

def getNumGradLJ(neighborList, latticeInfo, gradStep): # Calculate the gradient for LJ potential by numerical differentiation    
    neighborList_old = deepcopy(neighborList) # Make a copy of the original neighborList
    U_tot_LJ_old, U_avg_LJ_old = get_U_LJ(neighborList_old) # Get the original total energy
    delta = gradStep/100.0*latticeInfo[0] # Step is in unit of percentage of lattice constant; this computes the step size for the numerical differentiation
    
    gradLJ = []
    for i in range(len(neighborList_old)):
        atomGradLJ = []
        for j in range(1, 4): # Loop the three dimensions
            neighborList = deepcopy(neighborList_old) # Start with original neighborList

            coorChange = [0.0, 0.0, 0.0]
            coorChange[j - 1] = delta # Change the corresponding coordinate (x, y or z) by the numerical diff. step size

            # Update the information in the neighborList due to the finite difference change 
            neighborList = updateAtomCoor(neighborList, neighborList[i][0][4], coorChange)
            neighborList = updateAtomDist(neighborList)

            U_tot_LJ, U_avg_LJ = get_U_LJ(neighborList) # Calculate the new total energy due to the change in that coordinate of that atom
            atomGradLJ.append((U_tot_LJ - U_tot_LJ_old)/delta) # Partial derivative
        
        gradLJ.append(atomGradLJ)

    return neighborList_old, gradLJ

def getDescentStepBySecantLJ(neighborList, latticeInfo, gradLJ, gradStep, secantStep): # Get the descent step (alpha) by secant method for the LJ potential
    neighborList_old = deepcopy(neighborList)

    # Refer to Felix's slides on Lab 4
    # Calculate the nominator for the secant method: sigma*(del[f(x)] dot -del[f(x)])
    nominator = 0.0
    for i in range(len(gradLJ)):
        for j in range(len(gradLJ[i])):
            nominator += gradLJ[i][j]*(-gradLJ[i][j])
    nominator *= secantStep/100.0 # secantStep is in the unit of percentage of gradient

    # Calculate the denominator for the secant method: (del[f(x + delta)] - del[f(x)]) dot -del[f(x)]
    for i in range(len(neighborList)):
        atomID = neighborList[i][0][4]
        coorChange = [secantStep/100.0*-j for j in gradLJ[i]]
        neighborList = updateAtomCoor(neighborList, atomID, coorChange)
    neighborList = updateAtomDist(neighborList)
    
    neighborList, gradLJ_new = getNumGradLJ(neighborList, latticeInfo, gradStep) # del[f(x + delta)]
    
    gradLJ_diff = []
    for i in range(len(gradLJ)):
        gradLJ_diff.append([gradLJ_new[i][j] - gradLJ[i][j] for j in range(len(gradLJ[i]))]) # del[f(x + delta)] - del[f(x)]
    
    denominator = 0.0
    for i in range(len(gradLJ_diff)):
        for j in range(len(gradLJ_diff[i])):
            denominator += gradLJ_diff[i][j]*(-gradLJ[i][j])
    
    alpha = -nominator/denominator # alpha

    return neighborList_old, alpha

def getNumGradTer(neighborList, latticeInfo, gradStep, fileName): # Calculate the gradient for Tersoff potential by numerical differentiation    
    neighborList_old = deepcopy(neighborList) # Make a copy of the original neighborList
    U_tot_Ter_old, U_avg_Ter_old = get_U_Ter(neighborList_old, fileName) # Get the original total energy
    delta = gradStep/100.0*latticeInfo[0] # Step is in unit of percentage of lattice constant; this computes the step size for the numerical differentiation
    
    gradTer = []
    for i in range(len(neighborList_old)):
        atomGradTer = []
        for j in range(1, 4): # Loop the three dimensions
            neighborList = deepcopy(neighborList_old) # Start with original neighborList

            coorChange = [0.0, 0.0, 0.0]
            coorChange[j - 1] = delta # Change the corresponding coordinate (x, y or z) by the numerical diff. step size

            # Update the information in the neighborList due to the finite difference change 
            neighborList = updateAtomCoor(neighborList, neighborList[i][0][4], coorChange)
            neighborList = updateAtomDist(neighborList)

            U_tot_Ter, U_avg_Ter = get_U_Ter(neighborList, fileName) # Calculate the new total energy due to the change in that coordinate of that atom
            atomGradTer.append((U_tot_Ter - U_tot_Ter_old)/delta) # Partial derivative
        
        gradTer.append(atomGradTer)

    return neighborList_old, gradTer

def getDescentStepBySecantTer(neighborList, latticeInfo, gradTer, gradStep, secantStep, fileName): # Get the descent step (alpha) by secant method for the Tersoff potential
    neighborList_old = deepcopy(neighborList)

    # Refer to Felix's slides on Lab 4
    # Calculate the nominator for the secant method: sigma*(del[f(x)] dot -del[f(x)])
    nominator = 0.0
    for i in range(len(gradTer)):
        for j in range(len(gradTer[i])):
            nominator += gradTer[i][j]*(-gradTer[i][j])
    nominator *= secantStep/100.0 # secantStep is in the unit of percentage of gradient

    # Calculate the denominator for the secant method: (del[f(x + delta)] - del[f(x)]) dot -del[f(x)]
    for i in range(len(neighborList)):
        atomID = neighborList[i][0][4]
        coorChange = [secantStep/100.0*-j for j in gradTer[i]]
        neighborList = updateAtomCoor(neighborList, atomID, coorChange)
    neighborList = updateAtomDist(neighborList)
    
    neighborList, gradTer_new = getNumGradTer(neighborList, latticeInfo, gradStep, fileName) # del[f(x + delta)]
    
    gradTer_diff = []
    for i in range(len(gradTer)):
        gradTer_diff.append([gradTer_new[i][j] - gradTer[i][j] for j in range(len(gradTer[i]))]) # del[f(x + delta)] - del[f(x)]
    
    denominator = 0.0
    for i in range(len(gradTer_diff)):
        for j in range(len(gradTer_diff[i])):
            denominator += gradTer_diff[i][j]*(-gradTer[i][j])
    
    alpha = -nominator/denominator # alpha

    return neighborList_old, alpha

def updateNeighborList(neighborList, grad, alpha): # Update the atom positions (therefore the neighbor lists) by the step (-alpha*gradient)
    for i in range(len(neighborList)): # For each atom
        atomID = neighborList[i][0][4] # ID of that atom
        coorChange = [-alpha*j for j in grad[i]] # Coordinate change to be added to that atom's position
        neighborList = updateAtomCoor(neighborList, atomID, coorChange) # Update coordinates
    neighborList = updateAtomDist(neighborList) # Update separation distances

    return neighborList

def getNegGrad(grad): # Multiply all the gradient components by -1.0 
    for i in range(len(grad)):
        for j in range(len(grad[i])):
            grad[i][j] *= -1.0

    return grad

def getGamma(g_new, g_old): # Calculate gamma (refer to Felix's slides on Lab 4)
    nominator = 0.0
    denominator = 0.0
    for i in range(len(g_new)):
        for j in range(len(g_new[i])):
            nominator += g_new[i][j]*(g_new[i][j] - g_old[i][j]) # g_new dot (g_new - g_old)
            denominator += g_old[i][j]**2 # g_old dot g_old

    gamma = nominator/denominator
    if gamma > 0:
        return gamma
    else:
        return 0

def updateH(g_new, gamma, h_old): # Update h (refer to Felix's slides on Lab 4)
    for i in range(len(g_new)):
        for j in range(len(g_new[i])):
            h_old[i][j] = g_new[i][j] + gamma*h_old[i][j] # h_new = g_new + gamma*h_old

    return h_old

def printStepToXYZ(fileName, neighborList, U_stable, U_perturbed, U_step, i): # Print .xyz file of the steps for animation in VMD
    openFile = open(fileName + "_steps.xyz", "a")
  
    if i != 0:
        openFile.write("\n")

    openFile.write(str(len(neighborList)))
    openFile.write("\nFrame: %g; Minimum total energy: %gerg; Total energy at start (perturbed): %gerg; Total energy now: %gerg; Discrepancy from stable energy: %g%%" % (i, U_stable, U_perturbed, U_step, abs((U_step - U_stable)/U_stable*100.0)))
    for j in range(len(neighborList)):
            openFile.write("\n" + neighborList[j][0][0] + "\t" + str(neighborList[j][0][1]) + "\t" + str(neighborList[j][0][2]) + "\t" + str(neighborList[j][0][3]))

###################################### Script ######################################
################################# For LJ potential #################################
# Prepare a 2*2*2 Fcc Xenon lattice, with half of the lattice length as distance cutoff (parameters hard-coded into the functions)
# Structure Fcc; lattice constant 6.13A; periodicities 2*2*2; element Xe; dist. cutoff 6.129A; epsilon 320e-16erg; sigma 3.98A

# u_ prefix distinguishes from the variables used in functions
u_fileName = pAp1_LJ.get_xyzFile() # Generate .xyz file from user inputs
#u_distCutoff = float(raw_input("\nDistance cutoff for neighbor list (a positive number, in Angstrom; should be shorter than half of the length of the lattice): ")) # The distance cutoff
u_distCutoff = 6.129
print "Distance cutoff: %g Angstrom" % u_distCutoff

u_outputFileName = u_fileName[:-4] + "_neighborList" + ".txt"
u_neighborList, u_latticeInfo = getNeighborListFile(u_fileName, u_distCutoff, u_outputFileName) # The output neighbor lists file is created in the same directory as this .py file

# Format neighborList to the correct coordinates
u_neighborList = formatNeighborList(u_neighborList, u_latticeInfo)

# Perturb the lattice
u_neighborList_perturbed, u_neighborList_stable = getPerturbedList(u_neighborList, u_latticeInfo, 20.0) # Each coordinate is perturbed by 20.0% of lattice constant

# Calculate energies
U_tot_LJ_stable, U_avg_LJ_stable = get_U_LJ(u_neighborList_stable) # For the Fcc structure
U_tot_LJ_perturbed, U_avg_LJ_perturbed = get_U_LJ(u_neighborList_perturbed) # For the perturbed structure

######################## Steepest descent for LJ potential #########################
# Uncomment to use code
'''
discrepancy = 999.0 # discrepancy = |[(total energy now) - (minimum energy)]/(minimum energy)|
u_neighborList_new = deepcopy(u_neighborList_perturbed)
frame = 0 # Track the number of iteration

while discrepancy > 1.0: # Iteration stops when discrepancy between energy now and minimum energy is <= 1%
    U_tot_LJ_step, U_avg_LJ_step = get_U_LJ(u_neighborList_new) # Energies for the structure now
    discrepancy = abs((U_tot_LJ_step - U_tot_LJ_stable)/U_tot_LJ_stable*100.0) # Update discrepancy
    #print U_tot_LJ_stable, U_tot_LJ_perturbed, U_tot_LJ_step, discrepancy # Monitor

    printStepToXYZ("LJ_SD_" + u_fileName[:-4], u_neighborList_new, U_tot_LJ_stable, U_tot_LJ_perturbed, U_tot_LJ_step, frame) # Printing to output file

    u_neighborList_old, u_gradLJ = getNumGradLJ(u_neighborList_new, u_latticeInfo, 1.0) # Calculate the gradient. 1.0% of lattice constant used as step size
    u_neighborList_old, u_alpha = getDescentStepBySecantLJ(u_neighborList_old, u_latticeInfo, u_gradLJ, 1.0, 1.0) # Calculate alpha. 1.0% of lattice constant used as step size for gradient; 1.0% of gradient used as step size for secant method
    u_neighborList_new = updateNeighborList(u_neighborList_old, u_gradLJ, u_alpha) # Update the atoms' coordinates
    
    frame += 1
'''

####################### Conjugate gradient for LJ potential ########################
# Uncomment to use code
'''
discrepancy = 999.0 # discrepancy = |[(total energy now) - (minimum energy)]/(minimum energy)|
u_neighborList_new = deepcopy(u_neighborList_perturbed)
frame = 0 # Track the number of iteration

while discrepancy > 1.0: # Iteration stops when discrepancy between energy now and minimum energy is <= 1%
    U_tot_LJ_step, U_avg_LJ_step = get_U_LJ(u_neighborList_new) # Energies for the structure now
    discrepancy = abs((U_tot_LJ_step - U_tot_LJ_stable)/U_tot_LJ_stable*100.0) # Update discrepancy
    #print U_tot_LJ_stable, U_tot_LJ_perturbed, U_tot_LJ_step, discrepancy # Monitor
    
    printStepToXYZ("LJ_CG_" + u_fileName[:-4], u_neighborList_new, U_tot_LJ_stable, U_tot_LJ_perturbed, U_tot_LJ_step, frame) # Printing to output file

    u_neighborList_old, u_gradLJ = getNumGradLJ(u_neighborList_new, u_latticeInfo, 1.0) # Calculate the gradient. 1.0% of lattice constant used as step size
    
    if frame == 0: # Initialize the first value of g and h
        u_g_new = getNegGrad(u_gradLJ)
        u_h_new = deepcopy(u_g_new)
    else: # Starting from the second iteration
        u_g_dummy = getNegGrad(u_gradLJ) # g_dummy is g_new, but the g_new now is the g from last iteration, so use a dummmy variable
        u_gamma = -getGamma(u_g_dummy, u_g_new) # Get gamma from g_new and g_old
        u_g_new = deepcopy(u_g_dummy)
        u_h_new = updateH(u_g_new, u_gamma, u_h_new) # h_new = g_new + gamma*h_old

        u_gradLJ = getNegGrad(u_h_new) # Multiply all entries in h_new by -1.0; this is how I programmed the function
        
    u_neighborList_old, u_alpha = getDescentStepBySecantLJ(u_neighborList_old, u_latticeInfo, u_gradLJ, 1.0, 1.0) # Calculate alpha. 1.0% of lattice constant used as step size for gradient; 1.0% of gradient used as step size for secant method
    u_neighborList_new = updateNeighborList(u_neighborList_old, u_gradLJ, u_alpha) # Update the atoms' coordinates

    frame += 1
'''

############################### For Tersoff potential ##############################
# Prepare a 1*1*1 Diamond Silicon lattice, with half of the lattice length as distance cutoff (parameters hard-coded into the functions)
# Structure Diamond; lattice constant 5.43A; periodicities 1*1*1; element Si; dist. cutoff 2.714A; constants for the Tersoff potential are listed above

# u_ prefix distinguishes from the variables used in functions
u_fileName = pAp1_Ter.get_xyzFile() # Generate .xyz file from user inputs
#u_distCutoff = float(raw_input("\nDistance cutoff for neighbor list (a positive number, in Angstrom; should be shorter than half of the length of the lattice): ")) # The distance cutoff
u_distCutoff = 2.714
print "Distance cutoff: %g Angstrom" % u_distCutoff

u_outputFileName = u_fileName[:-4] + "_neighborList" + ".txt"
u_neighborList, u_latticeInfo = getNeighborListFile(u_fileName, u_distCutoff, u_outputFileName) # The output neighbor lists file is created in the same directory as this .py file
u_outputFileName = u_fileName[:-4] + ".xyz" # For the .xyz atom list file later used in functions

# Format neighborList to the correct coordinates
u_neighborList = formatNeighborList(u_neighborList, u_latticeInfo)

# Perturb the lattice
u_neighborList_perturbed, u_neighborList_stable = getPerturbedList(u_neighborList, u_latticeInfo, 20.0) # Each coordinate is perturbed by 20.0% of lattice constant

# Calculate energies
U_tot_Ter_stable, U_avg_Ter_stable = get_U_Ter(u_neighborList_stable, u_outputFileName) # For the diamond structure
U_tot_Ter_perturbed, U_avg_Ter_perturbed = get_U_Ter(u_neighborList_perturbed, u_outputFileName) # For the perturbed structure

###################### Steepest descent for Tersoff potential ######################
# Uncomment to use code
'''
discrepancy = 999.0 # discrepancy = |[(total energy now) - (minimum energy)]/(minimum energy)|
u_neighborList_new = deepcopy(u_neighborList_perturbed)
frame = 0 # Track the number of iteration

while discrepancy > 1.0: # Iteration stops when discrepancy between energy now and minimum energy is <= 1%
    U_tot_Ter_step, U_avg_Ter_step = get_U_Ter(u_neighborList_new, u_outputFileName) # Energies for the structure now
    discrepancy = abs((U_tot_Ter_step - U_tot_Ter_stable)/U_tot_Ter_stable*100.0) # Update discrepancy
    #print U_tot_Ter_stable, U_tot_Ter_perturbed, U_tot_Ter_step, discrepancy # Monitor
    
    printStepToXYZ("Ter_SD_" + u_fileName[:-4], u_neighborList_new, U_tot_Ter_stable, U_tot_Ter_perturbed, U_tot_Ter_step, frame) # Printing to output file

    u_neighborList_old, u_gradTer = getNumGradTer(u_neighborList_new, u_latticeInfo, 1.0, u_outputFileName) # Calculate the gradient. 1.0% of lattice constant used as step size
    u_neighborList_old, u_alpha = getDescentStepBySecantTer(u_neighborList_old, u_latticeInfo, u_gradTer, 1.0, 1.0, u_outputFileName) # Calculate alpha. 1.0% of lattice constant used as step size for gradient; 1.0% of gradient used as step size for secant method
    u_neighborList_new = updateNeighborList(u_neighborList_old, u_gradTer, u_alpha) # Update the atoms' coordinates

    frame += 1
'''

#################### Conjugate gradient for Tersoff potential ######################
# Uncomment to use code
'''
discrepancy = 999.0 # discrepancy = |[(total energy now) - (minimum energy)]/(minimum energy)|
u_neighborList_new = deepcopy(u_neighborList_perturbed)
frame = 0 # Track the number of iteration

while discrepancy > 1.0: # Iteration stops when discrepancy between energy now and minimum energy is <= 1%
    U_tot_Ter_step, U_avg_Ter_step = get_U_Ter(u_neighborList_new, u_outputFileName) # Energies for the structure now
    discrepancy = abs((U_tot_Ter_step - U_tot_Ter_stable)/U_tot_Ter_stable*100.0) # Update discrepancy
    #print U_tot_Ter_stable, U_tot_Ter_perturbed, U_tot_Ter_step, discrepancy # Monitor
    
    printStepToXYZ("Ter_CG_" + u_fileName[:-4], u_neighborList_new, U_tot_Ter_stable, U_tot_Ter_perturbed, U_tot_Ter_step, frame) # Printing to output file

    u_neighborList_old, u_gradTer = getNumGradTer(u_neighborList_new, u_latticeInfo, 1.0, u_outputFileName) # Calculate the gradient. 1.0% of lattice constant used as step size
    
    if frame == 0: # Initialize the first value of g and h
        u_g_new = getNegGrad(u_gradTer)
        u_h_new = deepcopy(u_g_new)
    else: # Starting from the second iteration
        u_g_dummy = getNegGrad(u_gradTer) # g_dummy is g_new, but the g_new now is the g from last iteration, so use a dummmy variable
        u_gamma = getGamma(u_g_dummy, u_g_new) # Get gamma from g_new and g_old
        u_g_new = deepcopy(u_g_dummy)
        u_h_new = updateH(u_g_new, u_gamma, u_h_new) # h_new = g_new + gamma*h_old

        u_gradTer = getNegGrad(u_h_new) # Multiply all entries in h_new by -1.0; this is how I programmed the function
        
    u_neighborList_old, u_alpha = getDescentStepBySecantTer(u_neighborList_old, u_latticeInfo, u_gradTer, 1.0, 1.0, u_outputFileName) # Calculate alpha. 1.0% of lattice constant used as step size for gradient; 1.0% of gradient used as step size for secant method
    u_neighborList_new = updateNeighborList(u_neighborList_old, u_gradTer, u_alpha) # Update the atoms' coordinates

    frame += 1
'''