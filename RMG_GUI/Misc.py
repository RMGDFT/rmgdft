# written by Wenchang Lu at NCSU


from PyQt5 import QtGui, QtCore, QtWidgets
from PyQt5 import QtWidgets as myQtW

from distutils.sysconfig import get_python_lib

## presumably need to import pickle module twice to make it work properly:
#import pickle
#import pickle
import json
import codecs

class Misc(QtWidgets.QWidget):
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
        myQtW.QWidget.__init__(self, parent)

        ################ __init__ : define the non-GUI variables: 

        ##################  __init__ : set up GUI:
        
        # Main layout
        self._layout = myQtW.QVBoxLayout()
        self.setLayout(self._layout)

        # Setup Groupbox
        group_box = myQtW.QGroupBox('Misc')
        self._layout.addWidget(group_box)

        form_layout = myQtW.QFormLayout()
        group_box.setLayout(form_layout)


        # Kohn Sham mgrid levels

        label = myQtW.QLabel('Kohn Sham MG level ')
        self._khlevel = myQtW.QSpinBox()
        self._khlevel.setValue(3)
        form_layout.addRow(label, self._khlevel)

        # poisson mgrid levels

        label = myQtW.QLabel('Poisson MG level ')
        self._poissonlevel = myQtW.QSpinBox()
        self._poissonlevel.setValue(4)
        form_layout.addRow(label, self._poissonlevel)


        # Kohn Sham time step

        label = myQtW.QLabel('Kohn Sham time step')
        self._khstep = myQtW.QLineEdit()
        validator = QtGui.QDoubleValidator(self._khstep)
        self._khstep.setValidator(validator)
        #form_layout.addRow(label, self._khstep)
        self._khstep.setText('0.66')

        # Poisson time step

        label = myQtW.QLabel('Poisson time step')
        self._poissonstep = myQtW.QLineEdit()
        validator = QtGui.QDoubleValidator(self._poissonstep)
        self._poissonstep.setValidator(validator)
        #form_layout.addRow(label, self._poissonstep)
        self._poissonstep.setText('0.66')

        # Kohn-sham Solver
        label = myQtW.QLabel('Kohn-Sham Solver')
        self.KS_solver = myQtW.QComboBox()
        self.KS_solver.addItems([
                            "davidson",
                            "multigrid"
                            ])
        form_layout.addRow(label, self.KS_solver)

        # VH Solver
        label = myQtW.QLabel('Poisson Solver')
        self.VH_solver = myQtW.QComboBox()
        self.VH_solver.addItems([
                            "pfft",
                            "multigrid"
                            ])
        form_layout.addRow(label, self.VH_solver)

        # Subdiag Solver
        label = myQtW.QLabel('Subdiag Solver')
        self.subdiag_solver = myQtW.QComboBox()
        self.subdiag_solver.addItems([
                            "cusolver",
                            "lapack",
                            "scalapack"
                            ])
        form_layout.addRow(label, self.subdiag_solver)

        self.local_pp_delocalization = myQtW.QCheckBox('local pseudopotentials extended to whole cell(delocalization)')
        form_layout.addRow(self.local_pp_delocalization)
        self.proj_delocalization = myQtW.QCheckBox('non-local projectors extended to whole cell(delocalization)')
        form_layout.addRow(self.proj_delocalization)
        self.folded_spectrum = myQtW.QCheckBox('use folded spectrum method')
        form_layout.addRow(self.folded_spectrum)

        label = myQtW.QLabel('states count and occupation ')
        self._state_count = myQtW.QLineEdit()
        form_layout.addRow(label, self._state_count)
        self._state_count.setText('')

        label = myQtW.QLabel('states count and occupation spin_up ')
        self._state_count_up = myQtW.QLineEdit()
        form_layout.addRow(label, self._state_count_up)
        self._state_count_up.setText('')

        label = myQtW.QLabel('states count and occupation spin_down')
        self._state_count_down = myQtW.QLineEdit()
        form_layout.addRow(label, self._state_count_down)
        self._state_count_down.setText('')

        label = myQtW.QLabel('potential acceleration constant')
        self._p_acc_const= myQtW.QLineEdit()
        validator = QtGui.QDoubleValidator(self._p_acc_const)
        self._p_acc_const.setValidator(validator)
        form_layout.addRow(label, self._p_acc_const)
        self._p_acc_const.setText('1.0')

        label = myQtW.QLabel('number of unoccupied states')
        self._num_unocc = myQtW.QSpinBox()
        self._num_unocc.setValue(10)
        form_layout.addRow(label, self._num_unocc)
    def state(self):
        """
           @return A dictionary containing the widget state.
        """
        input_misc_lines  =   '\n# ****  Multigrid **** \n\n'
        input_misc_lines +=   'kohn_sham_mg_levels  = "'               + str(self._khlevel.text()) +'"\n'
        input_misc_lines +=   'poisson_mg_levels    = "'               + str(self._poissonlevel.text()) +'"\n'
        #input_misc_lines +=   'kohn_sham_time_step  = "'               + str(self._khstep.text()) +'"\n'
        #input_misc_lines +=   'poisson_time_step    = "'               + str(self._poissonstep.text()) +'"\n'

        input_misc_lines += """

# RMG supports a pure multigrid Kohn-Sham solver as well as
# a multigrid preconditioned davidson solver. The davidson
# solver is usually better for smaller problems with the pure
# multigrid solver often being a better choice for very large
# problems.
#kohn_sham_solver="davidson"
#kohn_sham_solver="multigrid"

#poisson solver can be either pfft or multigrid
"""

        input_misc_lines +=   'kohn_sham_solver  = "'               + (self.KS_solver.currentText()) +'"\n'
        input_misc_lines +=   'poisson_solver  = "'               + (self.VH_solver.currentText()) +'"\n'

        input_misc_lines += """

# RMG supports a variety of subspace diagonalization options depending
# on the hardware and libraries available for a specific platform
"""

        input_misc_lines +=   'subdiag_diver  = "'               + (self.subdiag_solver.currentText()) +'"\n'

        if self._state_count.text() != '': 
            input_misc_lines +=   ' states_count_and_occupation = "'  + (self._state_count.text()) +'"\n'
        if self._state_count_up.text() != '': 
            input_misc_lines +=   ' states_count_and_occupation_spin_up = "'  + (self._state_count_up.text()) +'"\n'
        if self._state_count_down.text() != '': 
            input_misc_lines +=   ' states_count_and_occupation_spin_down = "'  + (self._state_count_down.text()) +'"\n'

        input_misc_lines += 'unoccupied_states_per_kpoint = "%s"'%self._num_unocc.text()


        input_misc_lines += """

# The Beta function projectors for a particular ion decay rapidly
# in real-space with increasing r. For large cells truncating the
# real-space representation of the projector can lead to
# significant computational savings with a small loss of accuracy.
# For smaller cells the computational cost is the same for localized
# and delocalized projectors so it is better to set localize_projectors
# to false.
#
# similar for local pseudopotential
"""

        if(self.proj_delocalization.isChecked()):
            input_misc_lines += 'localize_projectors = "false"\n'
        else:
            input_misc_lines += 'localize_projectors = "true"\n'

        if(self.local_pp_delocalization.isChecked()):
            input_misc_lines += 'localize_localpp = "false"\n'
        else:
            input_misc_lines += 'localize_localpp = "true"\n'
            
        if(self.folded_spectrum.isChecked()):
            input_misc_lines += 'folded_spectrum = "true"\n'

        input_misc_lines += 'potential_acceleration_constant_step = "'+ self._p_acc_const.text() + '"\n'


        state={ 'input_misc_lines': input_misc_lines }
        return state

    ######### end of state(self):


