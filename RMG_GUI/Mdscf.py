# written by Wenchang Lu at NCSU


from PyQt4 import QtGui, QtCore

from distutils.sysconfig import get_python_lib

## presumably need to import pickle module twice to make it work properly:
#import pickle
#import pickle
import json
import codecs

class Mdscf(QtGui.QWidget):
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
        group_box = QtGui.QGroupBox('SCF')
        self._layout.addWidget(group_box)

        form_layout = QtGui.QFormLayout()
        group_box.setLayout(form_layout)

        # max SCF steps

        label = QtGui.QLabel('Max scf steps for RMG ')
        self._maxscf = QtGui.QLineEdit()
        self._maxscf.setText('50')
        form_layout.addRow(label, self._maxscf)

        # RMS convergence critierion
        label = QtGui.QLabel('RMS Convergence Criterion')
        self._rms = QtGui.QLineEdit()
        validator = QtGui.QDoubleValidator(self._rms)
        self._rms.setValidator(validator)
        form_layout.addRow(label, self._rms)
        self._rms.setText('1e-10')

        label = QtGui.QLabel('energy Convergence Criterion')
        self._rms_energy = QtGui.QLineEdit()
        validator = QtGui.QDoubleValidator(self._rms)
        self._rms_energy.setValidator(validator)
        form_layout.addRow(label, self._rms_energy)
        self._rms_energy.setText('1e-12')

#        label = QtGui.QLabel('Max scf steps for NEGF ')
#        self._maxscf_NEGF= QtGui.QLineEdit()
#        self._maxscf_NEGF.setText('30')
#        form_layout.addRow(label, self._maxscf_NEGF)

        # RMS convergence critierion
#        label = QtGui.QLabel('RMS Convergence Criterion for NEGF')
#        self._rms_NEGF = QtGui.QLineEdit()
#        validator = QtGui.QDoubleValidator(self._rms_NEGF)
#        self._rms_NEGF.setValidator(validator)
#        form_layout.addRow(label, self._rms_NEGF)
#        self._rms_NEGF.setText('1e-12')


        # Setup Groupbox
        group_box = QtGui.QGroupBox('Relax/MD')
        self._layout.addWidget(group_box)

        form_layout = QtGui.QFormLayout()
        group_box.setLayout(form_layout)
        # max MD steps

        label = QtGui.QLabel('Max MD/FastRelax steps ')
        self._maxmd = QtGui.QLineEdit()
        self._maxmd.setText('10')
        form_layout.addRow(label, self._maxmd)

        # force convergence critierion
        label = QtGui.QLabel('Max force (Ha/au) < ')
        self._forcemax = QtGui.QLineEdit()
        validator = QtGui.QDoubleValidator(self._forcemax)
        self._forcemax.setValidator(validator)
        form_layout.addRow(label, self._forcemax)
        self._forcemax.setText('1e-3')

        # ionic time step

        label = QtGui.QLabel('ionic time step')
        self._ionstep = QtGui.QLineEdit()
        validator = QtGui.QDoubleValidator(self._ionstep)
        self._ionstep.setValidator(validator)
        form_layout.addRow(label, self._ionstep)
        self._ionstep.setText('45')

        # relax_dynamic_timestep
        label = QtGui.QLabel('time step changed dynamicall')
        self._rdt = QtGui.QComboBox()
        self._rdt.addItems(["true", "false"])
        form_layout.addRow(label, self._rdt)
    
        # Setup Groupbox
        group_box = QtGui.QGroupBox('Mixing')
        self._layout.addWidget(group_box)

        form_layout = QtGui.QFormLayout()
        group_box.setLayout(form_layout)

        # charge density mixing parameter

        label = QtGui.QLabel('Charge Density mixing')
        self._qmix = QtGui.QLineEdit()
        validator = QtGui.QDoubleValidator(self._qmix)
        self._qmix.setValidator(validator)
        form_layout.addRow(label, self._qmix)
        self._qmix.setText('0.1')

        # projector mixing parameter

        label = QtGui.QLabel('Projector mixing')
        self._pmix = QtGui.QLineEdit()
        validator = QtGui.QDoubleValidator(self._pmix)
        self._pmix.setValidator(validator)
        #form_layout.addRow(label, self._pmix)
        self._pmix.setText('0.1')

        #  charge density mixing method
        label = QtGui.QLabel('charge density mixing method')
        self._mixmethod = QtGui.QComboBox()
        self._mixmethod.addItems(["Pulay", "Linear", "Broyden"])
        form_layout.addRow(label, self._mixmethod)



        Hlayout = QtGui.QHBoxLayout()
        #Pulay mixing histroy
        label = QtGui.QLabel('  Pulay Order ')
        self._pulayorder = QtGui.QSpinBox()
        self._pulayorder.setValue(5)
        Hlayout.addWidget(label)
        Hlayout.addWidget(self._pulayorder)

        #Pulay mixing scale
        label = QtGui.QLabel(' scale (beta)  ')
        self._pulaybeta = QtGui.QLineEdit()
        self._pulaybeta.setText('0.5')
        Hlayout.addWidget(label)
        Hlayout.addWidget(self._pulaybeta)

        #Pulay mixing refresh
        label = QtGui.QLabel(' refresh steps  ')
        self._pulayrefresh = QtGui.QLineEdit()
        self._pulayrefresh.setText('100')
        Hlayout.addWidget(label)
        Hlayout.addWidget(self._pulayrefresh)
        form_layout.addRow(Hlayout)

        #K Point
        

    def state(self):
        """
           @return A dictionary containing the widget state.
        """

        input_mdscf_lines =   '\n# **** MD/Relax controls  **** \n  \n'
        input_mdscf_lines +=   'max_md_steps = "'               + str(self._maxmd.text()) +'"\n'
        input_mdscf_lines +=   'fast_relax_max_force = "'           + str(self._forcemax.text()) +'"\n'
        input_mdscf_lines +=   'ionic_time_step = "'           + str(self._ionstep.text()) +'"\n'
        input_mdscf_lines +=   'relax_dynamic_timestep = "'           + self._rdt.currentText() +'"\n'

        input_mdscf_lines +=   '\n# **** Mixing controls **** \n  \n'
        input_mdscf_lines +=   """

# RMG supports Broyden, Pulay and Linear mixing
# When the davidson Kohn-Sham solver is selected Broyden or
# Pulay are preferred. For the multigrid solver Linear with
# potential acceleration is often (but not always) the best
# choice.
#charge_mixing_type = "Broyden"
#charge_mixing_type = "Pulay"
#charge_mixing_type = "Linear"
"""

        input_mdscf_lines +=   'charge_mixing_type = "'           + self._mixmethod.currentText() +'"\n'
        input_mdscf_lines +=   'charge_density_mixing = "'      + str(self._qmix.text()) +'"\n'
        input_mdscf_lines +=   'charge_pulay_order = "'           + self._pulayorder.text() +'"\n'
        input_mdscf_lines +=   'charge_pulay_scale = "'           + self._pulaybeta.text() +'"\n'
        input_mdscf_lines +=   'charge_pulay_refresh = "'           + self._pulayrefresh.text() +'"\n'

        input_mdscf_lines +=   '\n# **** SCF controls ****  \n  \n'
        input_mdscf_lines +=   'max_scf_steps = "'              + str(self._maxscf.text()) +'"\n'
        input_mdscf_lines +=   'rms_convergence_criterion = "'           + str(self._rms.text()) +'"\n'
        input_mdscf_lines +=   'energy_convergence_criterion = "'           + str(self._rms_energy.text()) +'"\n'
#        input_mdscf_lines +=   input_mdscf_lines

#        input_mdscf_lines +=   '\n# **** SCF controls ****  \n  \n'
#        input_mdscf_lines +=   'max_scf_steps = "'              + str(self._maxscf_NEGF.text()) +'"\n'
#        input_mdscf_lines +=   'rms_convergence_criterion = "'           + str(self._rms_NEGF.text()) +'"\n'
#        input_mdscf_lines +=   input_mdscf_lines
        state={ 'input_mdscf_lines': input_mdscf_lines }
#                'input_mdscf_lines': input_mdscf_lines_NEGF}
        return state

    ######### end of state(self):


