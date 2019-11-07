# The Chinese University of Hong Kong, Department of Physics
# PHYS4061 Project A Part I
# Written by Yip Hoi Tung (1155092406)
# September 2018

# Get inputs from the user: which crystal structure, the lattice constant, the periodicities for all 3 directions
def get_inputs(structureList): # structureList is a list of the names of the crystal structures supported by this program
    structureCheck = False
    while structureCheck == False: # Check if the user's input for the crystal structure is in the supported list, otherwise prompt user to re-enter
        #structure = raw_input("\nCrystal structure (%s): " % "/".join(structureList))
        structure = "Diamond"
        print "Crystal structure: %s" % structure

        structureCheck = structure in structureList
        if structureCheck == False:
            print "Requested structure not found! Please input again.\n"

    #latticeC = float(raw_input("Lattice constant (a positive number, in Angstrom): ")) # Lattice constant
    latticeC = 5.43
    print "Lattice constant: %g Angstrom" % latticeC

    #periodX = int(raw_input("Periodicity in the x-direction (a positive integer): "))
    periodX = 1
    #periodY = int(raw_input("Periodicity in the y-direction (a positive integer): "))
    periodY = 1
    #periodZ = int(raw_input("Periodicity in the z-direction (a positive integer): "))
    periodZ = 1
    print "Periodicity: %g*%g*%g" % (periodX, periodY, periodZ)

    period = [periodX, periodY, periodZ] # A list of the periodicities

    #element = raw_input("Which element are the atoms? ") # Which element the atoms are
    element = "Si"
    print "Element: %s" % element

    return structure, latticeC, period, element

# Return the unit cell for that crystal structure
# The unit cell is a list of atoms; each atom is a list including its x, y and z-coordinate in the format: [element, x, y, z]
def get_cell(structure, latticeC, element):
    if structure == "Sc":
        return [[element, 0.0 * latticeC, 0.0 * latticeC, 0.0 * latticeC]]
    if structure == "Bcc":
        return [[element, 0.0 * latticeC, 0.0 * latticeC, 0.0 * latticeC], [element, 0.5 * latticeC, 0.5 * latticeC, 0.5 * latticeC]]
    if structure == "Fcc":
        return [[element, 0.0 * latticeC, 0.0 * latticeC, 0.0 * latticeC], [element, 0.0 * latticeC, 0.5 * latticeC, 0.5 * latticeC], [element, 0.5 * latticeC, 0.0 * latticeC, 0.5 * latticeC], [element, 0.5 * latticeC, 0.5 * latticeC, 0.0 * latticeC]]
    if structure == "Diamond":
        return [[element, 0.0 * latticeC, 0.0 * latticeC, 0.0 * latticeC], [element, 0.0 * latticeC, 0.5 * latticeC, 0.5 * latticeC], [element, 0.5 * latticeC, 0.0 * latticeC, 0.5 * latticeC], [element, 0.5 * latticeC, 0.5 * latticeC, 0.0 * latticeC], [element, 0.25 * latticeC, 0.25 * latticeC, 0.25 * latticeC], [element, 0.75 * latticeC, 0.25 * latticeC, 0.75 * latticeC], [element, 0.25 * latticeC, 0.75 * latticeC, 0.75 * latticeC], [element, 0.75 * latticeC, 0.75 * latticeC, 0.25 * latticeC]]
    if structure == "BC8":
        x_ref = 0.100299699 # Positional parameter for BC8 silicon
        x = float(raw_input("Input the positional parameter x for the BC8 structure (suggested value is %g): " % x_ref))
        return [[element, 0.0 * latticeC, 0.0 * latticeC, 0.0 * latticeC], [element, (1.0 - 2.0*x) * latticeC, (0.5 - 2.0*x) * latticeC, 0.0 * latticeC], [element, 0.0 * latticeC, (1.0 - 2.0*x) * latticeC, (0.5 - 2.0*x) * latticeC], [element, (0.5 - 2.0*x) * latticeC, 0.0 * latticeC, (1.0 - 2.0*x) * latticeC], [element, 0.5 * latticeC, 0.5 * latticeC, 0.5 * latticeC], [element, (0.5- 2.0*x) * latticeC, (0.5 - 2.0*x) * latticeC, (0.5 - 2.0*x) * latticeC], [element, 0.5 * latticeC, (0.5 - 2.0*x) * latticeC, (1.0 - 2.0*x) * latticeC], [element, (0.5 - 2.0*x) * latticeC, 0.5 * latticeC, 0.0 * latticeC], [element, (1.0 - 2.0*x) * latticeC, 0.5 * latticeC, (0.5 - 2.0*x) * latticeC], [element, 0.0 * latticeC, (0.5 - 2.0*x) * latticeC, 0.5 * latticeC], [element, (0.5 - 2.0*x) * latticeC, (1.0 - 2.0*x) * latticeC, 0.5 * latticeC], [element, 0.5 * latticeC, 0.0 * latticeC, (0.5 - 2.0*x) * latticeC], [element, 0.5 * latticeC, (1.0 - 2.0*x) * latticeC, 0.0 * latticeC], [element, (1.0 - 2.0*x) * latticeC, 0.0 * latticeC, 0.5 * latticeC], [element, 0.0 * latticeC, 0.5 * latticeC, (1.0 - 2.0*x) * latticeC], [element, (1.0 - 2.0*x) * latticeC, (1.0 - 2.0*x) * latticeC, (1.0 - 2.0*x) * latticeC]]

# Return the components for the XYZ file format: number of atoms, comment and a list of the atoms
def get_outputs(structure, latticeC, period, element):
    comment = "Crystal structure: %s; Lattice constant: %g Angstrom; Periodicity (x * y * z): %g * %g * %g" % (structure, latticeC, period[0], period[1], period[2])

    # The atom list is produced by repeating the unit cell in the x-direction, then repeating the resulted 'line' of atoms in the y-direction, then finally repeating the resulted 'slice' in the z-direction
    cell = get_cell(structure, latticeC, element)
    atomList = [x[:] for x in cell] # At first the atom list has only the atoms from the unit cell
    for d in range(3): # Loop for the 3 dimensions
        atomBlock = [x[:] for x in atomList] # Copy the 'block' to be repeated (either a unit cell, a 'line' or a 'slice' of atoms depending on the dimension being worked on)
        for i in range(1, period[d]): # The number of times of repeating by the periodicity
            for n in range(len(atomBlock)): # A repeated 'block' is created by moving each atom by the lattice constant in that direction
                atomBlock[n][d + 1] += latticeC
            atomList.extend([x[:] for x in atomBlock]) # Add the created block to the atom list

    # Check if the resulted list of atoms contains the correct number of atoms
    numOfAtoms = len(atomList)
    if numOfAtoms != len(cell) * period[0] * period[1] * period[2]:
        print 1/0

    return numOfAtoms, comment, atomList

def get_xyzFile():
    # The u_ prefix distinguishes themselves from the variables used in the functions
    u_structureList = ["Sc", "Bcc", "Fcc", "Diamond", "BC8"] # The list of crystal structures supported
    u_structure, u_latticeC, u_period, u_element = get_inputs(u_structureList)
    u_numOfAtoms, u_comment, u_atomList = get_outputs(u_structure, u_latticeC, u_period, u_element)

    # Producing the output file
    nameOfFile = "%s_%s_%g_%g_%g.xyz" % (u_element, u_structure, u_period[0], u_period[1], u_period[2]) # .xyz format
    outputFile = open(nameOfFile, "w")
    outputFile.write(str(u_numOfAtoms)) # First line is the number of atoms
    outputFile.write("\n" + u_comment) # Second line is the comment
    for i in range(len(u_atomList)): # The remaining lines are the atoms
        outputFile.write("\n" + "\t".join([str(x) for x in u_atomList[i]]))
    outputFile.close()
    print "The atom list (%s) has been generated." % nameOfFile

    return nameOfFile