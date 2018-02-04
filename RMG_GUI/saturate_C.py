#  calculate how many neighbors for each C atom and add hydrogen atoms
#  as needed. it works for graphene nanoribbon structures and not tested
#  for other systems
# written by Wenchang Lu at NCSU


from math import *
def saturate_C(elements, coords, a,b,c):

    CHbond = 1.101/0.529177
    bondcri = 1.80/0.529177
    neighbor=[]
    element_add = []
    coords_add = []
    for i in range(len(elements)):
        if(elements[i] == 'C'):
            v0 = coords[i]
            k = 0
            del neighbor
            neighbor=[]
            for j in range(len(elements)):
                v1 = coords[j]
                distance = norm_period(v0, v1, a,b,c)
                if((distance < bondcri) & (distance >0.0001) ): 
                    neighbor.append(j)
                    k+=1
            if( k > 4 or k <2 ): 
                print 'warning:  Carbon atom %d has %d neighbors'%(i,k)
                print v0
            if( k == 3): 
                v1 = coords[neighbor[0]]
                v2 = coords[neighbor[1]]
                v3 = coords[neighbor[2]]
                check_and_add3(v0,v1,v2,v3,element_add, coords_add, CHbond,a,b,c)
            if( k == 2): 
                v1 = coords[neighbor[0]]
                v2 = coords[neighbor[1]]
                check_and_add2(v0,v1,v2,element_add, coords_add, CHbond,a,b,c)
        if(elements[i] == 'H'):
            v0 = coords[i]
            k = 0
            for j in range(len(elements)):
                v1 = coords[j]
                distance = norm_period(v0, v1, a,b,c)
                if((distance < CHbond + 0.3) & (distance >0.0001) ): 
                    k+=1
            if (k !=1): 
                print 'warning:  Hydrogen atom %d has %d neighbors'%(i,k)
                
    for e in element_add: elements.append(e)
    for c in coords_add: coords.append(c)
    return coords_add

                
def check_and_add2(v0,v1,v2,element_add, coords_add, CHbond,a,b,c):
    vector1 = [v1[0] - v0[0],v1[1] - v0[1], v1[2] - v0[2]]
    vector2 = [v2[0] - v0[0],v2[1] - v0[1], v2[2] - v0[2]]
    unit_cell_fold(vector1, a,b,c)
    unit_cell_fold(vector2, a,b,c)
    norm(vector1)
    norm(vector2)
    angle = angle_2vectors(vector1, vector2)
    if(angle > 116.0):
        newvec = [vector1[0] + vector2[0],vector1[1] + vector2[1],vector1[2] + vector2[2]]
        norm(newvec)
        x = v0[0] - newvec[0] * CHbond
        y = v0[1] - newvec[1] * CHbond
        z = v0[2] - newvec[2] * CHbond
        newcoords = [x,y,z]
        element_add.append('H')
        coords_add.append(newcoords)
    else:
        print 'warning sp3 C not programed yet angle = %f'%angle
        
    return            

def check_and_add3(v0,v1,v2, v3,element_add, coords_add, CHbond, a,b,c):
    vector1 = [v1[0] - v0[0],v1[1] - v0[1], v1[2] - v0[2]]
    vector2 = [v2[0] - v0[0],v2[1] - v0[1], v2[2] - v0[2]]
    unit_cell_fold(vector1, a,b,c)
    unit_cell_fold(vector2, a,b,c)
    norm(vector1)
    norm(vector2)
    angle = angle_2vectors(vector1, vector2)
    
    if(angle <115): print 'warning sp3 C not programed yet angle = %f'%angle
    return            
    
def angle_2vectors(vec1, vec2):
    dis3 = vec1[0] * vec2[0] +vec1[1] * vec2[1] +vec1[2] * vec2[2] 
    angle = acos(dis3)/3.1415926 *180.0
    return angle

def norm_period(v0, v1, a,b,c):
    vector = [v0[0] - v1[0],v0[1] - v1[1], v0[2] - v1[2]]
    unit_cell_fold(vector, a,b,c)

    dis = vector[0] * vector[0] +vector[1] * vector[1] +vector[2] * vector[2] 
    dis = sqrt(dis)
    return dis
def unit_cell_fold(vector, a,b,c):
    if vector[0] > a/2.0: vector[0] -= a
    if vector[0] < -a/2.0: vector[0] += a

    if vector[1] > b/2.0: vector[1] -= b
    if vector[1] < -b/2.0: vector[1] += b
    
    if vector[2] > c/2.0: vector[2] -= c
    if vector[2] < -c/2.0: vector[2] += c
    return
def norm(vector):
    dis = vector[0] * vector[0] +vector[1] * vector[1] +vector[2] * vector[2] 
    dis = sqrt(dis)
    vector[0] = vector[0]/dis
    vector[1] = vector[1]/dis
    vector[2] = vector[2]/dis
    return
