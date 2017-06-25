#from numpy import *
import numpy  
import re
import CifFile

import pymol
from pymol import cmd
from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import QObject, pyqtSignal

from distutils.sysconfig import get_python_lib
from drawbox import *

## presumably need to import pickle module twice to make it work properly:
#import pickle
#import pickle
import json
import codecs

class conf:
    def __init__(self, lattice, elements, coords):
        self.lattice = lattice
        self.elements = elements
        self.coords = coords
class Configuration(QtGui.QWidget):
    """
       Widget for the basic setup of a NCSURMG calculation.
    """

    CONF_CHANGED = pyqtSignal()

    # member functions:

    def __init__(self, parent = None):
        """
           Constructor.

           @param parent : The parent widget.
        """

        ################ __init__ : Initialize base class 
        QtGui.QWidget.__init__(self, parent)

        ################ __init__ : define the non-GUI variables: 

        self.input_para={}
        
        elements = ['C']
        a = 12.234
        b = 14.234
        c = 8.156
        lattice = [a, b, c]

        coords = [ [ 4.078     ,  4.078     ,  0.        ],
                   [ 4.078     ,  8.156     ,  0.        ],
                   [ 4.078     ,  8.156     ,  0.        ],
                   [ 4.078     ,  8.156     ,  0.        ]]



        self.conf = []
        conf1 = conf(lattice, elements, coords)
        self.conf.append(conf1)
        self.conf.append(conf1)
        self.conf.append(conf1)


   # Main layout
        self._layout = QtGui.QVBoxLayout()
        self.setLayout(self._layout)


        group_box = QtGui.QGroupBox('   Choose input files for configuration   ')
        self._layout.addWidget(group_box)

        form_layout = QtGui.QFormLayout()
        group_box.setLayout(form_layout)

        self.layout  = QtGui.QHBoxLayout()

        self.label =   QtGui.QLabel('file for configuration')
        #self.inputfile =    QtGui.QLineEdit('/home/can/AMS-DATA.cif')
        #self.inputfile =    QtGui.QLineEdit('/home/can/projects/atk-on/RMG-GUI_readcif/input')
        self.inputfile =    QtGui.QLineEdit('./input')

        self.button =  QtGui.QPushButton('Browse...')
        self.pymolshow =  QtGui.QPushButton('Pymol Show')
        self.quitshow =  QtGui.QPushButton('reinit Pymol')

        self.layout.addWidget(self.label )
        self.layout.addWidget(self.inputfile  )
        self.layout.addWidget(self.button)
        self.layout.addWidget(self.pymolshow)
        self.layout.addWidget(self.quitshow)

        self.inputfile.setToolTip('edit the input file name')
        self.button.setToolTip('select an nput file name')

        self.button.clicked.connect( self.createConnect() )
#       self.createConnect()
        self.open_and_read(self.inputfile.text())
#       self.selectfile()
        self.pymolshow.clicked.connect( self.createshow() )
        self.quitshow.clicked.connect( self.quitpymol )

        form_layout.addRow(self.layout)
        
        self.lattlayout  = QtGui.QHBoxLayout()

#        lattlabel =   QtGui.QLabel('              lattice constant:')
       # self.latticea =  QtGui.QLineEdit('0.0   0.0   0.0')
#        self.latticea =  QtGui.QLabel('0.0   0.0   0.0')

#        self.inputfile.textChanged.connect(self.open_and_read(self.inputfile.text()))





###############padding and bbx#############

        self.bbxlayouts  = QtGui.QHBoxLayout()

        bbxlabel =   QtGui.QLabel('              Bounding Box:   a = ')
        self.bbxa =  QtGui.QLineEdit('1.0')
        self.bbxb =  QtGui.QLineEdit('1.0')
        self.bbxc =  QtGui.QLineEdit('1.0')

        self.bbxlayouts.addWidget(bbxlabel )
        self.bbxlayouts.addWidget(self.bbxa)
        self.bbxlayouts.addWidget(QtGui.QLabel('   b = '))
        self.bbxlayouts.addWidget(self.bbxb)
        self.bbxlayouts.addWidget(QtGui.QLabel('   c = '))
        self.bbxlayouts.addWidget(self.bbxc)

        form_layout.addRow(self.bbxlayouts)
        self.lattlayouts  = QtGui.QHBoxLayout()

        lattlabel =   QtGui.QLabel('              lattice constant:   a = ')
        self.latticea =  QtGui.QLineEdit('10.0')
        self.latticeb =  QtGui.QLineEdit('10.0')
        self.latticec =  QtGui.QLineEdit('10.0')

        self.lattlayouts.addWidget(lattlabel )
        self.lattlayouts.addWidget(self.latticea)
        self.lattlayouts.addWidget(QtGui.QLabel('   b = '))
        self.lattlayouts.addWidget(self.latticeb)
        self.lattlayouts.addWidget(QtGui.QLabel('   c = '))
        self.lattlayouts.addWidget(self.latticec)

        form_layout.addRow(self.lattlayouts)

        self.paddinglayouts  = QtGui.QHBoxLayout()
        self.paddingx = QtGui.QCheckBox('padding x')
        self.paddingy = QtGui.QCheckBox('padding y')
        self.paddingz = QtGui.QCheckBox('padding z')

        self.xpad =  QtGui.QLineEdit('10.0')
        self.ypad =  QtGui.QLineEdit('10.0')
        self.zpad =  QtGui.QLineEdit('10.0')

        self.paddinglayouts.addWidget(self.paddingx)
        self.paddinglayouts.addWidget(self.xpad)
        self.paddinglayouts.addWidget(self.paddingy)
        self.paddinglayouts.addWidget(self.ypad)
        self.paddinglayouts.addWidget(self.paddingz)
        self.paddinglayouts.addWidget(self.zpad)
        self.changepadding = QtGui.QPushButton('change padding')
        self.paddinglayouts.addWidget(self.changepadding)
        self.changepadding.clicked.connect( self.updatepadding)

        form_layout.addRow(self.paddinglayouts)

        
        self.latticea.textChanged.connect(self.updateunitcell)
        self.latticeb.textChanged.connect(self.updateunitcell)
        self.latticec.textChanged.connect(self.updateunitcell)
           
#        self.open_and_read()
#        self.inputfiles.textChanged.connect(self.open_and_read)

    def updateunitcell(self):
        self.conf.lattice[0]  = float(self.latticea.text())
        self.conf.lattice[1]  = float(self.latticeb.text())
        self.conf.lattice[2]  = float(self.latticec.text())
        self.CONF_CHANGED.emit()
        
    def updatepadding(self):
        ax = self.ax ;
        ay = self.ay ;
        az = self.az ;
        if(self.paddingx.isChecked()):
            ax = self.ax + float(self.xpad.text())
        if(self.paddingy.isChecked()):
            ay = self.ay + float(self.ypad.text())
        if(self.paddingz.isChecked()):
            az = self.az + float(self.zpad.text())
        self.latticea.setText(str(ax))
        self.latticeb.setText(str(ay))
        self.latticec.setText(str(az))

 

#        self.lattlayout.addWidget(lattlabel )
#        self.lattlayout.addWidget(self.latticea)


#        form_layout.addRow(self.lattlayout)
 

###############padding and bbx#############
    def selectfile(self):        

#        dialog = QtGui.QFileDialog(self)
#       dialog.exec_()
#        self.open_and_read(self.inputfile.text())  


        dialog = QtGui.QFileDialog(self)
        directory=os.getcwd()#QtGui.QFileDialog.getExistingDirectory(self)
        dialog.setDirectory(directory)# self.inputfile.text() )
        if QtCore.QDir( self.inputfile.text() ).exists():
            dialog.selectFile( self.inputfile.text() )
        if dialog.exec_():
            self.inputfile.setText(dialog.selectedFiles()[0])
        self.open_and_read(self.inputfile.text())  

    def createConnect(self):
        """
                   internal lambda-creating function to connect
           parameter-free signal from radio buttons for PP variants to
                   a parameter-accepting class method
        """
        return lambda: self.selectfile()


    def open_and_read(self,filename):
       filenames=filename.split('.')
       #filenames=re.split('.',filename)
       #print filenames[len(filenames)-2]
       #print filenames
       if filenames[len(filenames)-1]=='cif':
           self.open_and_read_cif(filename)
       if filenames[len(filenames)-1]=='xyz':
           self.open_and_read_xyz(filename)
#       if filenames[len(filenames)-1]=='input':
#           self.open_and_read_input(filename)
       if 'input' == filename.split('/')[len(filename.split('/'))-1]:
           self.open_and_read_input(filename)
           

    def open_and_read_input(self,filename):
       f = open (filename)
       input_file=f.readlines()
       input_para={}
       uncomment_input=[x for x in input_file if x[0]!='#' ]
       input_stream=''.join(uncomment_input)
       input_elements=input_stream.split('"')
       i=0
       for x in input_elements:
          if i==0:
              key=x.translate(None,'\n= ')
              i=1
          else:
              value=x.translate(None,'\n')
              i=0
              input_para.update({key:value})
#       print input_para
       self.input_para=input_para
       atoms=input_para['atoms']
#       print atoms
       conf_sp=atoms.split()
#       print conf_sp
       
       elements=[]
       lattice=[input_para['a_length'],input_para['b_length'],input_para['c_length']]
       coords=[]
       
       for i in range(len(conf_sp)/5):
          elements.append(conf_sp[i*5])
#          coord=conf_sp[i*5+1]+' '+conf_sp[i*5+2]+' '+conf_sp[i*5+3]
          coord=[float(conf_sp[i*5+1]),float(conf_sp[i*5+2]),float(conf_sp[i*5+3])]
          coords.append(coord)
       self.conf = conf(lattice, elements, coords)            
#           print elements
          


    def open_and_read_cif(self,filename):

        cf= CifFile.ReadCif(str(filename))
        self.input_para={}
         

        filekeys=cf.keys()

        cb=cf['global']
        keys=cb.keys()


        a=eval(cb['_cell_length_a'])  /0.529177
        b=eval(cb['_cell_length_b'])  /0.529177
        c=eval(cb['_cell_length_c'])  /0.529177
        pi = 3.1415926535
        alpha=eval(cb['_cell_angle_alpha'])/180.*pi
        beta=eval(cb['_cell_angle_beta'])/180.*pi
        gamma=eval(cb['_cell_angle_gamma'])/180.*pi
       

        symop=cb['_space_group_symop_operation_xyz']

        sites_base=cb.GetLoop('_atom_site_label')

            
        elements_base=cb['_atom_site_label']
#           print 'good'
#           occupancy_base=cb['_atom_site_occupancy']
        x_base=cb['_atom_site_fract_x']
        y_base=cb['_atom_site_fract_y']
        z_base=cb['_atom_site_fract_z']
#           print x_base
           
# site_base['_atom_site_fract_x']:


        i=0
        new_sites=[]
        for site in sites_base:
            newsite=[]
            for sym in symop:
                mapping={'x':x_base[i],'y':y_base[i],'z':z_base[i]}
                symeq=sym
                for k,v in mapping.iteritems():
                    symeq=symeq.replace(k,v)
                newcoord=[]
                for eq in symeq.split(','):
                    value=eval(eq)
                    if value==0:
                        value=0
                    newcoord.append(str(value))
                newelement=[elements_base[i],newcoord]
                if newelement not in newsite:#!=newcoord:
                    newsite.append(newelement)                     
            i+=1

#               symmatrix=i
            new_sites+=newsite


        coords=[]
        elements=[]
        lattice=[a,b,c]

        for site in new_sites:
            element= ''.join([i for i in site[0] if not i.isdigit()])
            elements.append(element)
         #elements.append(site[0])
            coords.append([float(site[1][0])*a,float(site[1][1])*b,float(site[1][2])*c])


        self.conf = conf(lattice, elements, coords)            


    def open_and_read_xyz(self, filename):
#        = self.inputfiles.text()
        try:
            a_to_bohr = 1.0/0.529177
            self.CONF_CHANGED.emit()
            f = open (filename)
            num_atoms = int(f.readline())
            comments = f.readline()
            elements = []
            coords = []
            for k in range(num_atoms):
                line = f.readline()
                linesp = line.split()
                elements.append(linesp[0])
                coords.append([float(linesp[1])*a_to_bohr, float(linesp[2])*a_to_bohr, float(linesp[3])*a_to_bohr])

            minl =  numpy.amin(coords, axis=0)
            maxl =  numpy.amax(coords, axis=0)
            self.ax = maxl[0] - minl[0] 
            self.ay = maxl[1] - minl[1]
            self.az = maxl[2] - minl[2]
            self.bbxa.setText(str(self.ax))
            self.bbxb.setText(str(self.ay))
            self.bbxc.setText(str(self.az))

            cellinfo = f.readline()
            line = f.readline()
            if line == '':
                print "no lattice information, unit cell is padding with 10A"
                
                lattice = [self.ax + 10.0, self.ay + 10.0, self.az + 10.0]
                self.conf = conf(lattice, elements, coords)            
                self.latticea.setText(str(self.ax + 10.0))
                self.latticeb.setText(str(self.ay + 10.0))
                self.latticec.setText(str(self.az + 10.0))
                    
                
            else:
                linesp = line.split()
                lattice = [float(linesp[0])*a_to_bohr, float(linesp[1])*a_to_bohr, float(linesp[2])*a_to_bohr]
   
                self.ax = float(linesp[0])*a_to_bohr
                self.ay = float(linesp[1])*a_to_bohr
                self.az = float(linesp[2])*a_to_bohr

                self.conf = conf(lattice, elements, coords)            
                self.latticea.setText(str(self.ax))
                self.latticeb.setText(str(self.ay))
                self.latticec.setText(str(self.az))
		print lattice
        except:
            print "no such a file === ", filename

    def createshow(self):
        """
                   internal lambda-creating function to connect
           parameter-free signal from radio buttons for PP variants to
                   a parameter-accepting class method
        """
        return lambda: self.structureshow()
    def structureshow(self):
#        pymol.finish_launching()
        pymol.finish_launching()
        filename = self.inputfile.text()
        cmd.load(str(filename))
        cmd.show('spheres')
        minx = 0.0
        miny = 0.0
        minz = 0.0
        maxx = self.conf.lattice[0]
        maxy = self.conf.lattice[1]
        maxz = self.conf.lattice[2]
        
        cmd.extend("drawbox", drawbox(minx, miny, minz, maxx, maxy, maxz, 2.0, 1,1,1))

    def quitpymol(self):
        cmd.reinitialize()

