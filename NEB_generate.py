#!/usr/bin/python

import os
if __name__ == '__main__':

    num_images = 7
    num_proc = 2
    first_image = "image00"
    last_image = "image06"
    
    print "modify lines in NEB_generate.py to generate NEB inputs for different images"
    ctrl_init = """ 
spin_polarization="false"
# spin_polarization="true"

image_per_node="1"
max_neb_steps = "10"
neb_spring_constant = "0.50000000"

#  max number of images is defined in params.h MAX_IMGS 99
num_images="%d"

# for each image, path is either relative to the path in job file or
# full path
image_infos="
./%s/  input   %d 
"""%(num_images,first_image,num_proc)

    image_per_node="1"
    max_neb_steps = "10"
    neb_spring_constant = "0.50000000"

    first_coor = []
    f = open(first_image+"/input")
    all_lines = f.readlines()
    f.close()
    input_keep ="" 
    for i in range(0, len(all_lines)):
        line =  all_lines[i]
        input_keep += line
        if "atoms" in line and "=" in line and "\"" in line:
            line_atom = i
            break;
    for i in range(line_atom+1, len(all_lines)):
        coor = all_lines[i]
        if "\"" in coor: break
        first_coor.append(coor)

#    print first_coor
            
    last_coor = []
    f = open(last_image+"/input")
    all_lines = f.readlines()
    f.close()
    for i in range(0, len(all_lines)):
        line =  all_lines[i]
        if "atoms" in line and "=" in line and "\"" in line:
            line_atom = i
            break;
    for i in range(line_atom+1, len(all_lines)):
        coor = all_lines[i]
        if "\"" in coor: break
        last_coor.append(coor)

    if len(first_coor) != len(last_coor): 
        print "first and last image have different number of atoms"
        exit()
    
    for i in range(1, num_images-1):

        pecent = 1.0/(num_images -1) * i
        input_image = input_keep
        if i <10: dir_name = "image0" + str(i)
        else : dir_name = "image" + str(i)
        if not os.path.exists(dir_name): os.mkdir(dir_name)
        ctrl_init += "./%s/  input   %d\n"%(dir_name, num_proc) 
        
        for j in range(0, len(first_coor)):
            atom1 = first_coor[j].split()    
            atom2 = last_coor[j].split()    
        
            atom3 = atom1
            x = float(atom1[1]) * (1.0-pecent) + float(atom2[1]) * pecent
            y = float(atom1[2]) * (1.0-pecent) + float(atom2[2]) * pecent
            z = float(atom1[3]) * (1.0-pecent) + float(atom2[3]) * pecent
            #print x,y,z, pecent
            atom3[1] = str(x)
            atom3[2] = str(y)
            atom3[3] = str(z)
            input_image += "  ".join(atom3) +"\n"
            #print "  ".join(atom3) +"\n"
        input_image += "\"\n"

        filename = dir_name + "/input"
        with open(filename, "w") as f: f.write(input_image)

    ctrl_init += "./%s/  input   %d\n\"\n"%(last_image, num_proc) 
    with open("ctrl_init.dat", "w") as f: f.write(ctrl_init)

