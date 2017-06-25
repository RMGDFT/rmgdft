# written by Wenchang Lu at NCSU

import os

# here is the function that will return the lattice constant from bravias_lattice vector
def lattice(lattice):
    """
    @return ascii string with lattice information 
    """

    # put parameters part into the snippet
    #calculate length of unit cell vectors


    a = lattice[0]
    b = lattice[1]
    c = lattice[2]

    # put unitcell information into the snippet

    snippet = """

# Lattice constants in (x, y, z) directions
# a, b, c, cos(alpha), cos(beta), cos(gamma)
a_length="%s" b_length="%s" c_length="%s"
alpha="0.0" beta="0.0" gamma="0.0"

""" %(a,b,c)

    return snippet


