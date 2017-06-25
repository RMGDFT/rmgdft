

from PyQt4 import QtGui, QtCore

from distutils.sysconfig import get_python_lib

## presumably need to import pickle module twice to make it work properly:
#import pickle
#import pickle
import json
import string
import codecs

class Setup(QtGui.QWidget):
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
        group_box = QtGui.QGroupBox('Setup')
        self._layout.addWidget(group_box)

        form_layout = QtGui.QFormLayout()
        group_box.setLayout(form_layout)

        #description
        label = QtGui.QLabel('Description of input')
        self._description = QtGui.QLineEdit()
        form_layout.addRow(label, self._description)
        self._description.setText('Short description of the input file')

        # Start mode 
        label = QtGui.QLabel('Start mode')
        self._start_mode = QtGui.QComboBox()
        self._start_mode.addItems([
                            "LCAO Start",
                            "Restart From File",
                            "Random Start"
                            ])
        form_layout.addRow(label, self._start_mode)



        # Calculation mode 
        label = QtGui.QLabel('calculation mode')
        self._calculation_mode = QtGui.QComboBox()
        self._calculation_mode.addItems([
                            "Quench Electrons",
                            "Relax Structure",
                            "NEB Relax",
                            "Constant Volume And Energy",
                            "Constant Temperature And Energy",
#                            "Constant Pressure And Energy(not implemented)",
                            "Constrained Fast Relax",
                            "Band Structure Only"
                             ])
        form_layout.addRow(label, self._calculation_mode)

        # Relax Method
        label = QtGui.QLabel('Relax Method')
        self._relaxmethod = QtGui.QComboBox()
        self._relaxmethod.addItems(["Fast Relax",
                "FIRE", "Quick Min", "MD Min", "LBFGS"])
        form_layout.addRow(label, self._relaxmethod)


        # XC selector:

        label = QtGui.QLabel('Exchange correlation functional')
        self._ecf = QtGui.QComboBox()
        self._ecf.addItems(['AUTO_XC', 'LDA' , 'GGA BLYP', 'GGA XB CP', 'GGA PBE'
                            ])
        self._ecf.setToolTip(self.tr('Select the exchange-correlation functional.'))
        form_layout.addRow(label, self._ecf)

        #bravais_lattice_type 

        label = QtGui.QLabel('Bravais lattice')
        self._brav = QtGui.QComboBox()
        self._brav.addItems([
                            'Orthorhombic Primitive',
                            'Cubic Primitive',
                            'Hexagonal Primitive'
                            ])
        form_layout.addRow(label, self._brav)

        # Charge
        label = QtGui.QLabel('System Charge')
        self._charge = QtGui.QLineEdit()
        validator = QtGui.QDoubleValidator(self._charge)
        self._charge.setValidator(validator)
        form_layout.addRow(label, self._charge)
        self._charge.setText('0.0')


        # Spin polarized
        self._spin_polarized = QtGui.QCheckBox('Spin polarized')
        form_layout.addRow(self._spin_polarized)

        # units
        group_box = QtGui.QGroupBox('Units')
        self._layout.addWidget(group_box)

        form_layout = QtGui.QFormLayout()
        group_box.setLayout(form_layout)

        #Length units
        label = QtGui.QLabel('Length unit')
        self._lengthunit = QtGui.QComboBox()
        self._lengthunit.addItems(["Bohr","Angstrom"])
        form_layout.addRow(label, self._lengthunit)

        # Atomic coordinate type
        label = QtGui.QLabel('Atomic Coordinate')
        self._atomcoor = QtGui.QComboBox()
        self._atomcoor.addItems(["Absolute","Cell Relative"])
        form_layout.addRow(label, self._atomcoor)


        # Occupation 
        group_box = QtGui.QGroupBox('Occupation')
        self._layout.addWidget(group_box)

        form_layout = QtGui.QFormLayout()
        group_box.setLayout(form_layout)


        #Occupation type
        label = QtGui.QLabel('Occupation type')
        self._occ = QtGui.QComboBox()
        self._occ.addItems(['Fermi Dirac','Fixed' 
                            ])
        form_layout.addRow(label, self._occ)

        #  Occupation temperature 
        label = QtGui.QLabel('electron temperature (eV)')
        self._occtem = QtGui.QLineEdit()
        validator = QtGui.QDoubleValidator(self._occtem)
        self._occtem.setValidator(validator)
        self._occtem.setText('0.1')

        form_layout.addRow(label, self._occtem)


        #  Occupation temperature 
        label = QtGui.QLabel('Occupation mixing')
        self._occmix = QtGui.QLineEdit()
        validator = QtGui.QDoubleValidator(self._occmix)
        self._occmix.setValidator(validator)
        self._occmix.setText('0.3')
        form_layout.addRow(label, self._occmix)


        # Machines
        group_box = QtGui.QGroupBox('Computer and project name')
        self._layout.addWidget(group_box)

        form_layout = QtGui.QFormLayout()
        group_box.setLayout(form_layout)

        label = QtGui.QLabel('Machine')
        self._machine = QtGui.QComboBox()
        self._machine.addItems(['bluewater','ox','Jaguar','chugach' ])
        form_layout.addRow(label, self._machine)

#        label = QtGui.QLabel('Project name')
#        self._projname = QtGui.QComboBox()
#        self._projname.addItems(['CHM022','ONRDC17053C5A','ONRDC17053240'])
#        form_layout.addRow(label, self._projname)
	label=QtGui.QLabel('Project name')
	self._projname=QtGui.QLineEdit()
	validator=QtGui.QDoubleValidator(self._projname)
	self._projname.setValidator(validator)
	self._projname.setText('test')


        label = QtGui.QLabel('queue type')
        self._queue = QtGui.QComboBox()
        self._queue.addItems(['batch','debug' ])
        form_layout.addRow(label, self._queue)




    def state(self):
        """
           @return A dictionary containing the widget state.
        """
        input_setup_lines = '\n# **** Setup  ****\n\n'
        input_setup_lines += 'description = "'                + self._description.text() +'"\n'
        input_setup_lines += 'start_mode = "'                + self._start_mode.currentText() +'"\n'
        input_setup_lines += 'calculation_mode = \"' +self._calculation_mode.currentText() +'\"\n'
        input_setup_lines += 'relax_method =\"' +self._relaxmethod.currentText()+'\"\n'
        input_setup_lines += 'exchange_correlation_type = "' + self._ecf.currentText() +'"\n'
        input_setup_lines += 'bravais_lattice_type = "' + self._brav.currentText() +'"\n'
        input_setup_lines += 'system_charge = "'             + str(self._charge.text()) +'"\n'
        input_occupation_lines  = '\n# *****Occupation *****\n\n'
        input_occupation_lines += 'occupations_type = "'           + self._occ.currentText() +'"\n'
	input_occupation_lines += 'occupation_electron_temperature_eV = "'             + str(self._occtem.text()) +'"\n'
        input_occupation_lines += 'occupation_number_mixing = "'             + str(self._occmix.text()) +'"\n'

        input_units_lines = '\n# **** Units  ****\n\n'
        input_units_lines += 'length_units = "'                + self._lengthunit.currentText() +'"\n'
        input_units_lines += 'atomic_coordinate_type = "'                + self._atomcoor.currentText() +'"\n'

        state={ 'input_setup_lines': input_setup_lines,
                'input_units_lines': input_units_lines,
                'input_occupation_lines':input_occupation_lines,
                'length_units': str(self._lengthunit.currentText())
                 }
        return state

    ######### end of state(self):


