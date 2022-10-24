import numpy as np
import matplotlib.pyplot as plt

element_delta=[
['H', 0,0,'greenyellow'],
['He',17,0,'springgreen'],
['Li',0,1,'violet'],
['Be',1,1,'sandybrown'],
['B',12,1,'lavender'],
['C',13,1,'greenyellow'],
['N',14,1,'greenyellow'],
['O',15,1,'greenyellow'],
['F',16,1,'gold'],
['Ne',17,1,'springgreen'],
['Na',0,2,'violet'],
['Mg',1,2,'sandybrown'],
['Al',12,2,'dodgerblue'],
['Si',13,2,'lavender'],
['P',14,2,'greenyellow'],
['S',15,2,'greenyellow'],
['Cl',16,2,'gold'],
['Ar',17,2,'springgreen'],
['K' ,0,3,'violet'],
['Ca',1,3,'sandybrown'],
['Sc',2,3,'cyan'],
['Ti',3,3,'cyan'],
['V' ,4,3,'cyan'],
['Cr',5,3,'cyan'],
['Mn',6,3,'cyan'],
['Fe',7,3,'cyan'],
['Co',8,3,'cyan'],
['Ni',9,3,'cyan'],
['Cu',10,3,'cyan'],
['Zn',11,3,'cyan'],
['Ga',12,3,'dodgerblue'],
['Ge',13,3,'lavender'],
['As',14,3,'lavender'],
['Se',15,3,'greenyellow'],
['Br',16,3,'gold'],
['Kr',17,3,'springgreen'],
['Rb',0,4,'violet'],
['Sr',1,4,'sandybrown'],
['Y' ,2,4,'cyan'],
['Zr',3,4,'cyan'],
['Nb',4,4,'cyan'],
['Mo',5,4,'cyan'],
['Tc',6,4,'cyan'],
['Ru',7,4,'cyan'],
['Rh',8,4,'cyan'],
['Pd',9,4,'cyan'],
['Ag',10,4,'cyan'],
['Cd',11,4,'cyan'],
['In',12,4,'dodgerblue'],
['Sn',13,4,'dodgerblue'],
['Sb',14,4,'lavender'],
['Te',15,4,'lavender'],
['I' ,16,4,'gold'],
['Xe',17,4,'springgreen'],
['Cs',0,5,'violet'],
['Ba',1,5,'sandybrown'],
['Lu',2,5,'cyan'],
['Hf',3,5,'cyan'],
['Ta',4,5,'cyan'],
['W' ,5,5,'cyan'],
['Re',6,5,'cyan'],
['Os',7,5,'cyan'],
['Ir',8,5,'cyan'],
['Pt',9,5,'cyan'],
['Au',10,5,'cyan'],
['Hg',11,5,'cyan'],
['Tl',12,5,'dodgerblue'],  
['Pb',13,5,'dodgerblue'],  
['Bi',14,5,'dodgerblue'],
['Po',15,5,'lavender'],
['At',16,5,'gold'],
['Rn',17,5,'springgreen'],
]    

f=open('Delta-out.txt','r')
all_lines = f.readlines()
f.close()

delta={'At':'    ---    '}
for line in all_lines:
    if('#' in line): continue
    if('-' in line): continue
    linesp = line.split() 
    if('N' in linesp[1]): 
        delta.update({linesp[0]:'    ---    '})
    else:
        
        delta.update({linesp[0]:' '+linesp[1]+' '})

fig, ax = plt.subplots(figsize=(20,12),dpi=300)

d_h = 1.500

ax.axis("off")
# these are matplotlib.patch.Patch properties
#props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

# place a text box in upper left in axes coords
for ele in element_delta:
    x = ele[1] * 0.05
    y = 1.0-ele[2] * 0.10
    textstr = '\n'.join((
        r'',
        ele[0],
        r'%5s' %(delta[ele[0]],),
        r''))
    #textstr = '\n'.join((
    #    ele[0],
    #    r' %.3f ' % (d_h, ),
    #    r''))
    props = dict(facecolor=ele[3], alpha=0.0)
    ax.text(x, y, textstr, transform=ax.transAxes, fontsize=14,
        horizontalalignment="center",verticalalignment='top', bbox=props)

plt.show()
plt.savefig('delta.png')


