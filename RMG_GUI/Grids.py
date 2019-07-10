# written by Wenchang Lu at NCSU

from PyQt4 import QtGui, QtCore

from distutils.sysconfig import get_python_lib

from procgrid3d_to_2d import *

## presumably need to import pickle module twice to make it work properly:
import json
import codecs
from numpy import *

class Grids(QtGui.QWidget):
    """
       Widget for the setup grids, including processor grids and space grids
    """

    def __init__(self, parent = None):
        """
           Constructor.

           @param parent : The parent widget.
        """

        QtGui.QWidget.__init__(self, parent)

        try:
            # Main layout
            self._layout = QtGui.QVBoxLayout()
            self.setLayout(self._layout)
            
            # Setup Groupbox
            group_box = QtGui.QGroupBox('Grid information for RealSpace Multigrid')
            self._layout.addWidget(group_box)
            
            form_layout = QtGui.QFormLayout()
            group_box.setLayout(form_layout)

            layout = QtGui.QHBoxLayout()
            label = QtGui.QLabel('         choose approximate grid spacing:')
            self.gridspacing = QtGui.QLineEdit()
            validator = QtGui.QDoubleValidator(self.gridspacing)
            self.gridspacing.setValidator(validator)
            self.button = QtGui.QPushButton('update gridspacing')

            layout.addWidget(self.gridspacing)
            layout.addWidget(self.button)

            form_layout.addRow(label,layout)

            self.gridspacing.setText('0.35')

            Hlayout = QtGui.QHBoxLayout()

            Hlayout = QtGui.QHBoxLayout()
            label = QtGui.QLabel('   wave function grids:')
            self._Nx = QtGui.QSpinBox()
            self._Nx.setMaximum(9000)
            self._Nx.setValue(48)
            self._Ny = QtGui.QSpinBox()
            self._Ny.setMaximum(9000)
            self._Ny.setValue(48)
            self._Nz = QtGui.QSpinBox()
            self._Nz.setMaximum(9000)
            self._Nz.setValue(48)

            Hlayout.addWidget(QtGui.QLabel('     Nx='))
            Hlayout.addWidget(self._Nx)
            Hlayout.addWidget(QtGui.QLabel('     Ny='))
            Hlayout.addWidget(self._Ny)
            Hlayout.addWidget(QtGui.QLabel('     Nz='))
            Hlayout.addWidget(self._Nz)

            form_layout.addRow(label,Hlayout)

            label = QtGui.QLabel('   exact grid spacing:')
            self.exacthxyz = QtGui.QLabel('      hx = 0.35             hy = 0.35             hz = 0.35 bohr')

            form_layout.addRow(label,self.exacthxyz)


            label = QtGui.QLabel('   grid spacing anisotropy:')
            self.anisotropy = QtGui.QLabel('      0.0%')

            form_layout.addRow(label,self.anisotropy)

            Hlayout = QtGui.QHBoxLayout()
            label = QtGui.QLabel('   pot-rho/wave grid ratio:')
            self._potratio = QtGui.QSpinBox()
            self._potratio.setMaximum(8)
            self._potratio.setMinimum(1)
            self._potratio.setValue(2)
            Hlayout.addWidget(self._potratio)

            form_layout.addRow(label,Hlayout)

            Hlayout = QtGui.QHBoxLayout()
            label = QtGui.QLabel('   processor grids for 3D space:')
            self._Pex = QtGui.QSpinBox()
            self._Pex.setMaximum(9000)
            self._Pex.setValue(1)
            self._Pey = QtGui.QSpinBox()
            self._Pey.setMaximum(9000)
            self._Pey.setValue(1)
            self._Pez = QtGui.QSpinBox()
            self._Pez.setMaximum(9000)
            self._Pez.setValue(1)

            Hlayout.addWidget(QtGui.QLabel('    Pex='))
            Hlayout.addWidget(self._Pex)
            Hlayout.addWidget(QtGui.QLabel('    Pey='))
            Hlayout.addWidget(self._Pey)
            Hlayout.addWidget(QtGui.QLabel('    Pez='))
            Hlayout.addWidget(self._Pez)

            form_layout.addRow(label,Hlayout)

            # Setup Groupbox
            group_box = QtGui.QGroupBox('K-point ')
            self._layout.addWidget(group_box)

            form_layout = QtGui.QFormLayout()
            group_box.setLayout(form_layout)

            Hlayout = QtGui.QHBoxLayout()
            label = QtGui.QLabel('   k-point mesh ')
            self._kx = QtGui.QSpinBox()
            self._kx.setMaximum(9000)
            self._kx.setValue(1)
            self._ky = QtGui.QSpinBox()
            self._ky.setMaximum(9000)
            self._ky.setValue(1)
            self._kz = QtGui.QSpinBox()
            self._kz.setMaximum(9000)
            self._kz.setValue(1)
            self._kx.setMinimum(1)
            self._ky.setMinimum(1)
            self._kz.setMinimum(1)

            Hlayout.addWidget(QtGui.QLabel('     Kx='))
            Hlayout.addWidget(self._kx)
            Hlayout.addWidget(QtGui.QLabel('     Ky='))
            Hlayout.addWidget(self._ky)
            Hlayout.addWidget(QtGui.QLabel('     Kz='))
            Hlayout.addWidget(self._kz)

            form_layout.addRow(label,Hlayout)

            Hlayout = QtGui.QHBoxLayout()
            label = QtGui.QLabel('   k-point shift ')
            self._is_shift_x = QtGui.QSpinBox()
            self._is_shift_x.setMaximum(1)
            self._is_shift_x.setValue(0)
            self._is_shift_y = QtGui.QSpinBox()
            self._is_shift_y.setMaximum(1)
            self._is_shift_y.setValue(0)
            self._is_shift_z = QtGui.QSpinBox()
            self._is_shift_z.setMaximum(1)
            self._is_shift_z.setValue(0)

            Hlayout.addWidget(QtGui.QLabel('     is_shift_x='))
            Hlayout.addWidget(self._is_shift_x)
            Hlayout.addWidget(QtGui.QLabel('     is_shift_y='))
            Hlayout.addWidget(self._is_shift_y)
            Hlayout.addWidget(QtGui.QLabel('     is_shift_z='))
            Hlayout.addWidget(self._is_shift_z)

            form_layout.addRow(label,Hlayout)

            group_box = QtGui.QGroupBox('K-point for Band structure')
            self._layout.addWidget(group_box)

            form_layout = QtGui.QFormLayout()
            group_box.setLayout(form_layout)

            Hlayout = QtGui.QHBoxLayout()
            label = QtGui.QLabel('   number of special lines:')
            self.num_klines = QtGui.QSpinBox()
            self.num_klines.setValue(1)
            Hlayout.addWidget(self.num_klines)
            form_layout.addRow(label,Hlayout)
            label = QtGui.QLabel('   kx (2pi/a)             ky (2pi/b)                     kz (2pi/c)    number of kpoint    symbol:')
            form_layout.addRow(label)
            max_klines = 5
            self.kk_layouts = range(max_klines+1)
            self.kx_band = range(max_klines+1)
            self.ky_band = range(max_klines+1)
            self.kz_band = range(max_klines+1)
            self.kpts_band = range(max_klines+1)
            self.ksymbol_band = range(max_klines+1)
            for i in range(max_klines):
                self.kk_layouts[i]  = QtGui.QHBoxLayout() 
                self.kx_band[i] =  QtGui.QLineEdit() 
                self.kx_band[i].setText('0.0')
                self.kk_layouts[i].addWidget(self.kx_band[i] )
                self.ky_band[i] =  QtGui.QLineEdit() 
                self.ky_band[i].setText('0.0')
                self.kk_layouts[i].addWidget(self.ky_band[i] )
                self.kz_band[i] =  QtGui.QLineEdit() 
                self.kz_band[i].setText('0.0')
                self.kk_layouts[i].addWidget(self.kz_band[i] )
                self.kpts_band[i] =  QtGui.QSpinBox() 
                self.kpts_band[i].setValue(10)
                self.kk_layouts[i].addWidget(self.kpts_band[i] )
                self.ksymbol_band[i] =  QtGui.QLineEdit() 
                self.ksymbol_band[i].setText('G')
                self.kk_layouts[i].addWidget(self.ksymbol_band[i] )
                form_layout.addRow(self.kk_layouts[i])

        except:
            print " Grids layout error"

        try:
           self.button.clicked.connect(self.changeNxgrid)
           self._Nx.valueChanged.connect(self.changeothergrid)
           self._Ny.valueChanged.connect(self.changeNyzgrid)
           self._Nz.valueChanged.connect(self.changeNyzgrid)

           self.get2dgrid()

           self._Pex.valueChanged.connect(self.get2dgrid)
           self._Pey.valueChanged.connect(self.get2dgrid)
           self._Pez.valueChanged.connect(self.get2dgrid)

        except:
           print " Grids value change error"
    def changeNxgrid(self):
        Nxgrid = int(self.a/float(self.gridspacing.text()))
        Nxgrid = Nxgrid/self.gridfactor * self.gridfactor
        self._Nx.setValue(Nxgrid)

    def get2dgrid(self):
        
        pex = self._Pex.value()
        pey = self._Pey.value()
        pez = self._Pez.value()
        nprow, npcol = procgrid3d_to_2d(pex, pey, pez)
        self._nprow = nprow
        self._npcol = npcol
        

    def lattparameters(self, configuration, misc):

        vectors1 = configuration.conf.lattice
        self.a = float(vectors1[0] )
        self.b = float(vectors1[1])
        self.c = float(vectors1[2])
        khlevel = int(misc._khlevel.text())
        self.gridfactor = 2**khlevel
#        self.gridfactor = 1
#        for i in range(khlevel):
#            self.gridfactor *=2
        Nxgrid = int(self.a/float(self.gridspacing.text()))
        Nxgrid = Nxgrid/self.gridfactor * self.gridfactor
        self._Nx.setValue(Nxgrid)
        self._Nx.setSingleStep(self.gridfactor)
        self._Ny.setSingleStep(self.gridfactor)
        self._Nz.setSingleStep(self.gridfactor)
        self.changeothergrid()


    def changeNyzgrid(self):
#        try:
            
            hx = self.a/self._Nx.value()
            hy = self.b/self._Ny.value()
            hz = self.c/self._Nz.value()
            hmax = hx
            if hmax < hy: hmax = hy
            if hmax < hz: hmax = hz
            hmin = hx
            if hmin > hy: hmin = hy
            if hmin > hz: hmin = hz
    
            
            anis = str("%.1f" %((hmax/hmin-1) * 100))+'%'
            self.anisotropy.setText('      ' + anis )

            hx = str("%.4f" %hx)
            hy = str("%.4f" %hy)
            hz = str("%.4f" %hz)

            self.exacthxyz.setText('     hx = '+hx+'          hy = '+hy+'          hz = '+hz+' bohr')
    def changeothergrid(self):
#        try:
            
            ratio = self._Nx.value() /self.a 
            Nygrid = int(ratio * self.b)
            Nygrid = Nygrid/self.gridfactor * self.gridfactor
            item = (self.b -Nygrid /ratio ) * ratio
            if(item > self.gridfactor/2): Nygrid += self.gridfactor
            self._Ny.setValue(Nygrid)
            Nzgrid = int(ratio * self.c)
            Nzgrid = Nzgrid/self.gridfactor * self.gridfactor
            item = (self.c -Nzgrid /ratio ) * ratio
            if(item > self.gridfactor/2): Nzgrid += self.gridfactor
            self._Nz.setValue(Nzgrid)

            hx = self.a/self._Nx.value()
            hy = self.b/self._Ny.value()
            hz = self.c/self._Nz.value()
            hmax = hx
            if hmax < hy: hmax = hy
            if hmax < hz: hmax = hz
            hmin = hx
            if hmin > hy: hmin = hy
            if hmin > hz: hmin = hz
    
            
            anis = str("%.1f" %((hmax/hmin-1) * 100))+'%'
            self.anisotropy.setText('      ' + anis )

            hx = str("%.4f" %hx)
            hy = str("%.4f" %hy)
            hz = str("%.4f" %hz)

            self.exacthxyz.setText('     hx = '+hx+'          hy = '+hy+'          hz = '+hz+' bohr')



#        except:
#            print "load the coordinate files .xyz ... first"

    def state(self):
        """
           @return A dictionary containing the widget state.
        """
        try:
           input_grids_lines = '\n#wavefunction grid  and processor grid\n'
           input_grids_lines += 'wavefunction_grid ="' 
           input_grids_lines +='%d ' %self._Nx.value()
           input_grids_lines +='%d ' %self._Ny.value()
           input_grids_lines +='%d' %self._Nz.value() + '"\n'
           num_proc = self._Pex.value() * self._Pey.value() * self._Pez.value() 
           self.num_proc = num_proc
        #   input_grids_lines +='num_processor ="' +'%d'%num_proc +'"\n'

           input_grids_lines +='processor_grid="'
           input_grids_lines +='%d ' %self._Pex.value()
           input_grids_lines +='%d ' %self._Pey.value()
           input_grids_lines +='%d' %self._Pez.value() + '"\n'

           input_grids_lines +="""

# Ratio of the potential grid density to the wavefunction grid
# density. For example if the wavefunction grid is (72,72,72) and 
# potential_grid_refinement = "2" then the potential grid would be
# (144,144,144). The default value is 2 but it may sometimes be
# beneficial to adjust this. (For USPP the minimum value is also 2
# and it cannot be set lower. NCPP can be set to 1).
"""

           input_grids_lines +='potential_grid_refinement="'
           input_grids_lines +='%d ' %self._potratio.value() +'"\n'

           input_grids_lines +="""

#kpoint mesh set up
#kpoint_is_shift ="0 0 0" include gamma point
"""

           input_grids_lines += 'kpoint_mesh ="' 
           input_grids_lines +='%d ' %self._kx.value()
           input_grids_lines +='%d ' %self._ky.value()
           input_grids_lines +='%d' %self._kz.value() + '"\n'

           input_grids_lines += 'kpoint_is_shift ="' 
           input_grids_lines +='%d ' %self._is_shift_x.value()
           input_grids_lines +='%d ' %self._is_shift_y.value()
           input_grids_lines +='%d' %self._is_shift_z.value() + '"\n'

           input_grids_lines += 'kpoints_bandstructure ="\n' 
           for i in range(self.num_klines.value()+1):
               kx = self.kx_band[i].text()
               ky = self.ky_band[i].text()
               kz = self.kz_band[i].text()
               num = self.kpts_band[i].value()
               symbol = self.ksymbol_band[i].text()

               input_grids_lines += '%s  %s  %s  %d  %s\n'%(kx,ky,kz,num, symbol)
           input_grids_lines += '"\n'

    
        #   input_grids_lines +='Hamiltonia_processor_grid ="'
        #   input_grids_lines +='%d ' %self._nprow
        #   input_grids_lines +='%d' %self._npcol + '"\n'


        except:
           print "Grid  state error1"
           
        state={ 
                'input_grids_lines': input_grids_lines,
                }
        return state

    ######### end of state(self):


