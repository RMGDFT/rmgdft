# written by Wenchang Lu at NCSU


from PyQt4 import QtGui, QtCore

from distutils.sysconfig import get_python_lib

## presumably need to import pickle module twice to make it work properly:
#import pickle
#import pickle
import json
import codecs

class Misc(QtGui.QWidget):
    """
       Widget for the basic setup of a NCSURMG calculation.
    """


    # location of the plugin parts:
    _ncsurmg_addon_path = "/AddOns/NCSURMG/"


    # member functions:

    def __init__(self, parent = None):
        """
           Constructor.

           @param parent : The parent widget.
        """

        ################ __init__ : Initialize base class 
        QtGui.QWidget.__init__(self, parent)

        ################ __init__ : define the non-GUI variables: 

        ##################  __init__ : set up GUI:
        
        # Main layout
        self._layout = QtGui.QVBoxLayout()
        self.setLayout(self._layout)

        # Setup Groupbox
        group_box = QtGui.QGroupBox('Misc')
        self._layout.addWidget(group_box)

        form_layout = QtGui.QFormLayout()
        group_box.setLayout(form_layout)


        # Kohn Sham mgrid levels

        label = QtGui.QLabel('Kohn Sham MG level ')
        self._khlevel = QtGui.QSpinBox()
        self._khlevel.setValue(2)
        form_layout.addRow(label, self._khlevel)

        # poisson mgrid levels

        label = QtGui.QLabel('Poisson MG level ')
        self._poissonlevel = QtGui.QSpinBox()
        self._poissonlevel.setValue(2)
        form_layout.addRow(label, self._poissonlevel)


        # Kohn Sham time step

        label = QtGui.QLabel('Kohn Sham time step')
        self._khstep = QtGui.QLineEdit()
        validator = QtGui.QDoubleValidator(self._khstep)
        self._khstep.setValidator(validator)
        form_layout.addRow(label, self._khstep)
        self._khstep.setText('0.66')

        # Poisson time step

        label = QtGui.QLabel('Poisson time step')
        self._poissonstep = QtGui.QLineEdit()
        validator = QtGui.QDoubleValidator(self._poissonstep)
        self._poissonstep.setValidator(validator)
        form_layout.addRow(label, self._poissonstep)
        self._poissonstep.setText('0.66')

        # Kohn-sham Solver
        label = QtGui.QLabel('Kohn-Sham Solver')
        self.KS_solver = QtGui.QComboBox()
        self.KS_solver.addItems([
                            "multigrid",
                            "davidson"
                            ])
        form_layout.addRow(label, self.KS_solver)

        # VH Solver
        label = QtGui.QLabel('Poisson Solver')
        self.VH_solver = QtGui.QComboBox()
        self.VH_solver.addItems([
                            "pfft",
                            "multigrid"
                            ])
        form_layout.addRow(label, self.VH_solver)

        self.local_pp_delocalization = QtGui.QCheckBox('local pseudopotentials extended to whole cell(delocalization)')
        form_layout.addRow(self.local_pp_delocalization)
        self.proj_delocalization = QtGui.QCheckBox('non-local projectors extended to whole cell(delocalization)')
        form_layout.addRow(self.proj_delocalization)


    def state(self):
        """
           @return A dictionary containing the widget state.
        """
        input_misc_lines  =   '\n# ****  Multigrid **** \n\n'
        input_misc_lines +=   'kohn_sham_mg_levels  = "'               + str(self._khlevel.text()) +'"\n'
        input_misc_lines +=   'poisson_mg_levels    = "'               + str(self._poissonlevel.text()) +'"\n'
        input_misc_lines +=   'kohn_sham_time_step  = "'               + str(self._khstep.text()) +'"\n'
        input_misc_lines +=   'poisson_time_step    = "'               + str(self._poissonstep.text()) +'"\n'
        input_misc_lines +=   'kohn_sham_solver  = "'               + (self.KS_solver.currentText()) +'"\n'
        input_misc_lines +=   'poisson_solver  = "'               + (self.VH_solver.currentText()) +'"\n'

        if(self.local_pp_delocalization.isChecked()):
            input_misc_lines += 'localize_localpp = "false"\n'
        if(self.proj_delocalization.isChecked()):
            input_misc_lines += 'localize_projectors = "false"\n'


        state={ 'input_misc_lines': input_misc_lines }
        return state

    ######### end of state(self):


