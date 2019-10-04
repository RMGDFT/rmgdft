

from PyQt5 import QtGui, QtCore, QtWidgets
from PyQt5 import QtWidgets as myQtW

from distutils.sysconfig import get_python_lib

## presumably need to import pickle module twice to make it work properly:
#import pickle
#import pickle
import json
import string
import codecs

class Setup(myQtW.QWidget):
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
        group_box = myQtW.QGroupBox('Setup')
        self._layout.addWidget(group_box)

        form_layout = myQtW.QFormLayout()
        group_box.setLayout(form_layout)

        #description
        label = myQtW.QLabel('Description of input')
        self._description = myQtW.QLineEdit()
        form_layout.addRow(label, self._description)
        self._description.setText('Short description of the input file')

        # Start mode 
        label = myQtW.QLabel('Start mode')
        self._start_mode = myQtW.QComboBox()
        self._start_mode.addItems([
                            "LCAO Start",
                            "Restart From File",
                            "Random Start"
                            ])
        form_layout.addRow(label, self._start_mode)



        # Calculation mode 
        label = myQtW.QLabel('calculation mode')
        self._calculation_mode = myQtW.QComboBox()
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
        label = myQtW.QLabel('Relax Method')
        self._relaxmethod = myQtW.QComboBox()
        self._relaxmethod.addItems(["Fast Relax",
                "FIRE", "Quick Min", "MD Min", "LBFGS"])
        form_layout.addRow(label, self._relaxmethod)


        # XC selector:

        label = myQtW.QLabel('Exchange correlation functional')
        self._ecf = myQtW.QComboBox()
        self._ecf.addItems(['AUTO_XC', 'LDA' , 'GGA BLYP', 'GGA XB CP', 'GGA PBE'
                            ])
        self._ecf.setToolTip(self.tr('Select the exchange-correlation functional.'))
        form_layout.addRow(label, self._ecf)

        #bravais_lattice_type 

        label = myQtW.QLabel('Bravais lattice')
        self._brav = myQtW.QComboBox()
        self._brav.addItems([
                            'Orthorhombic Primitive',
                            'Cubic Primitive',
                            'Hexagonal Primitive'
                            ])
        form_layout.addRow(label, self._brav)

        # Charge
        label = myQtW.QLabel('System Charge')
        self._charge = myQtW.QLineEdit()
        validator = QtGui.QDoubleValidator(self._charge)
        self._charge.setValidator(validator)
        form_layout.addRow(label, self._charge)
        self._charge.setText('0.0')


        # Spin polarized
        self._spin_polarized = myQtW.QCheckBox('Spin polarized')
        form_layout.addRow(self._spin_polarized)

        # units
        group_box = myQtW.QGroupBox('Units')
        self._layout.addWidget(group_box)

        form_layout = myQtW.QFormLayout()
        group_box.setLayout(form_layout)

        #Length units
        label = myQtW.QLabel('Length unit')
        self._lengthunit = myQtW.QComboBox()
        self._lengthunit.addItems(["Bohr","Angstrom"])
        form_layout.addRow(label, self._lengthunit)

        # Atomic coordinate type
        label = myQtW.QLabel('Atomic Coordinate')
        self._atomcoor = myQtW.QComboBox()
        self._atomcoor.addItems(["Absolute","Cell Relative"])
        form_layout.addRow(label, self._atomcoor)


        # Occupation 
        group_box = myQtW.QGroupBox('Occupation')
        self._layout.addWidget(group_box)

        form_layout = myQtW.QFormLayout()
        group_box.setLayout(form_layout)


        #Occupation type
        label = myQtW.QLabel('Occupation type')
        self._occ = myQtW.QComboBox()
        self._occ.addItems(['Fermi Dirac','Fixed' 
                            ])
        form_layout.addRow(label, self._occ)

        #  Occupation temperature 
        label = myQtW.QLabel('electron temperature (eV)')
        self._occtem = myQtW.QLineEdit()
        validator = QtGui.QDoubleValidator(self._occtem)
        self._occtem.setValidator(validator)
        self._occtem.setText('0.1')

        form_layout.addRow(label, self._occtem)


        #  Occupation temperature 
        label = myQtW.QLabel('Occupation mixing')
        self._occmix = myQtW.QLineEdit()
        validator = QtGui.QDoubleValidator(self._occmix)
        self._occmix.setValidator(validator)
        self._occmix.setText('1.0')
        form_layout.addRow(label, self._occmix)


        # Machines
        group_box = myQtW.QGroupBox('Computer and project name')
        self._layout.addWidget(group_box)

        form_layout = myQtW.QFormLayout()
        group_box.setLayout(form_layout)

        label = myQtW.QLabel('Machine')
        self._machine = myQtW.QComboBox()
        self._machine.addItems(['Summit','bluewater'])
        form_layout.addRow(label, self._machine)
        label=myQtW.QLabel('Project name')

        self._projname=myQtW.QLineEdit()
        validator=QtGui.QDoubleValidator(self._projname)
        self._projname.setValidator(validator)
        self._projname.setText('CHP107')
        form_layout.addRow(label, self._projname)


        label = myQtW.QLabel('queue type')
        self._queue = myQtW.QComboBox()
        self._queue.addItems(['batch','debug' ])
        form_layout.addRow(label, self._queue)




    def state(self):
        """
           @return A dictionary containing the widget state.
        """
        input_setup_lines = '# description of run\n'
        input_setup_lines += 'description = "'                + self._description.text() +'"\n'

        input_setup_lines += """

# Uncommenting this will print out more information
#verbose="true"


# In most cases LCAO or Restart but certain special scenarios
# may require a Random or Modified LCAO start
#start_mode="LCAO Start"
#start_mode="Random Start"
#start_mode="Modified LCAO Start"
#start_mode="Restart From File"
"""
        input_setup_lines += 'start_mode = "'                + self._start_mode.currentText() +'"\n'

        input_setup_lines += """

# This is not an exhaustive list of options but does
# contain the most frequently used ones.
#calculation_mode="Quench Electrons"
#calculation_mode="Relax Structure"
#calculation_mode="Constant Volume And Energy"
#calculation_mode="Constant Temperature And Energy"
#calculation_mode="Band Structure Only"

"""

        input_setup_lines += 'calculation_mode = \"' +self._calculation_mode.currentText() +'\"\n'

        input_setup_lines += """

"""

        input_setup_lines += 'relax_method =\"' +self._relaxmethod.currentText()+'\"\n'

        input_setup_lines += """

# Most pseudopotentials specify the exchange correlation type they
# were generated with and the default value of AUTO_XC means that
# the type specified in the pseudopotial is what RMG will use. That
# can be overridden by specifying a value here. For a full list of
# the available types look in the source distribution at the file
# Headers/InputOpts.h around line 146.
"""

        input_setup_lines += 'exchange_correlation_type = "' + self._ecf.currentText() +'"\n'


        input_setup_lines += """

# RMG supports the following lattice types (Hexagonal at gamma-only)
#bravais_lattice_type="Cubic Primitive"
#bravais_lattice_type="Orthorhombic Primitive"
#bravais_lattice_type="Hexagonal Primitive"
"""

        input_setup_lines += 'bravais_lattice_type = "' + self._brav.currentText() +'"\n'
        input_setup_lines += 'system_charge = "'             + str(self._charge.text()) +'"\n'

        input_occupation_lines = """

# RMG supports several different ways of specifying orbital occupations.
# For a spin polarized system one may specify the occupations for up and
# down separately. In the case of a non-zero electronic temperature these
# will be adjusted as the calculation proceeds. For a non-spin polarized
# calculation look at one of the other examples.
#occupations_type = "Fixed"
#occupations_type = "Fermi Dirac"
#occupations_type = "MethfesselPaxton"
"""
        input_occupation_lines += 'occupations_type = "'           + self._occ.currentText() +'"\n'
        input_occupation_lines += 'occupation_electron_temperature_eV = "'             + str(self._occtem.text()) +'"\n'
        input_occupation_lines += 'occupation_number_mixing = "'             + str(self._occmix.text()) +'"\n'

        input_units_lines = """

#  length unit in Bohr or Angerstrom
# Default is Bohr
#crds_units="Bohr"
#crds_units="Angstrom"
"""

        input_units_lines += 'crds_units = "'                + self._lengthunit.currentText() +'"\n'


        input_units_lines += """

#atomic_coordinate_type="Cell Relative"
#atomic_coordinate_type="Absolute"
"""


        input_units_lines += 'atomic_coordinate_type = "'                + self._atomcoor.currentText() +'"\n'

        state={ 'input_setup_lines': input_setup_lines,
                'input_units_lines': input_units_lines,
                'input_occupation_lines':input_occupation_lines,
                'length_units': str(self._lengthunit.currentText())
                 }
        return state

    ######### end of state(self):



