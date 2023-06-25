from PyQt5 import QtGui
from PyQt5 import QtWidgets as myQtW

_NCSMURG_ADDON_PATH = "/AddOns/NCSURMG/"


class Setup(myQtW.QWidget):
    """
    Widget for the basic setup of a NCSURMG calculation.
    """

    # location of the plugin parts:
    _ncsurmg_addon_path = _NCSMURG_ADDON_PATH

    # member functions:

    def __init__(self, parent=None):
        """
        Constructor.

        @param parent : The parent widget.
        """

        # __init__ : Initialize base class
        myQtW.QWidget.__init__(self, parent)

        # __init__ : define the non-GUI variables:

        #  __init__ : set up GUI:

        # Main layout
        self._layout = myQtW.QVBoxLayout()
        self.setLayout(self._layout)

        # Setup Groupbox
        group_box = myQtW.QGroupBox("Setup")
        self._layout.addWidget(group_box)

        form_layout = myQtW.QFormLayout()
        group_box.setLayout(form_layout)

        # description
        label = myQtW.QLabel("Description of input")
        self._description = myQtW.QLineEdit()
        form_layout.addRow(label, self._description)
        self._description.setText("Short description of the input file")

        # Start mode
        label = myQtW.QLabel("Start mode")
        self._start_mode = myQtW.QComboBox()
        self._start_mode.addItems(
            [
                "LCAO Start",
                "Restart From File",
                "Random Start",
            ]
        )
        form_layout.addRow(label, self._start_mode)

        # Calculation mode
        label = myQtW.QLabel("calculation mode")
        self._calculation_mode = myQtW.QComboBox()
        self._calculation_mode.addItems(
            [
                "Quench Electrons",
                "Relax Structure",
                "NEB Relax",
                "Constant Volume And Energy",
                "Constant Temperature And Energy",
                "Constrained Fast Relax",
                "Band Structure Only",
            ]
        )
        form_layout.addRow(label, self._calculation_mode)

        # Relax Method
        label = myQtW.QLabel("Relax Method")
        self._relaxmethod = myQtW.QComboBox()
        self._relaxmethod.addItems(
            [
                "Fast Relax",
                "FIRE",
                "Quick Min",
                "MD Min",
                "LBFGS",
            ]
        )
        form_layout.addRow(label, self._relaxmethod)

        # XC selector:

        label = myQtW.QLabel("Exchange correlation functional")
        self._ecf = myQtW.QComboBox()
        self._ecf.addItems(
            [
                "AUTO_XC",
                "LDA",
                "GGA BLYP",
                "GGA XB CP",
                "GGA PBE",
            ]
        )
        self._ecf.setToolTip(
            self.tr("Select the exchange-correlation functional."),
        )
        form_layout.addRow(label, self._ecf)

        # bravais_lattice_type

        label = myQtW.QLabel("Bravais lattice")
        self._brav = myQtW.QComboBox()
        self._brav.addItems(
            [
                "Orthorhombic Primitive",
                "Cubic Primitive",
                "Hexagonal Primitive",
            ]
        )
        form_layout.addRow(label, self._brav)

        # Charge
        label = myQtW.QLabel("System Charge")
        self._charge = myQtW.QLineEdit()
        validator = QtGui.QDoubleValidator(self._charge)
        self._charge.setValidator(validator)
        form_layout.addRow(label, self._charge)
        self._charge.setText("0.0")

        # Spin polarized
        self._spin_polarized = myQtW.QCheckBox("Spin polarized")
        form_layout.addRow(self._spin_polarized)

        # units
        group_box = myQtW.QGroupBox("Units")
        self._layout.addWidget(group_box)

        form_layout = myQtW.QFormLayout()
        group_box.setLayout(form_layout)

        # Length units
        label = myQtW.QLabel("Length unit")
        self._lengthunit = myQtW.QComboBox()
        self._lengthunit.addItems(["Bohr", "Angstrom"])
        form_layout.addRow(label, self._lengthunit)

        # Atomic coordinate type
        label = myQtW.QLabel("Atomic Coordinate")
        self._atomcoor = myQtW.QComboBox()
        self._atomcoor.addItems(["Absolute", "Cell Relative"])
        form_layout.addRow(label, self._atomcoor)

        # Occupation
        group_box = myQtW.QGroupBox("Occupation")
        self._layout.addWidget(group_box)

        form_layout = myQtW.QFormLayout()
        group_box.setLayout(form_layout)

        # Occupation type
        label = myQtW.QLabel("Occupation type")
        self._occ = myQtW.QComboBox()
        self._occ.addItems(["Fermi Dirac", "Fixed"])
        form_layout.addRow(label, self._occ)

        #  Occupation temperature
        label = myQtW.QLabel("electron temperature (eV)")
        self._occtem = myQtW.QLineEdit()
        validator = QtGui.QDoubleValidator(self._occtem)
        self._occtem.setValidator(validator)
        self._occtem.setText("0.1")

        form_layout.addRow(label, self._occtem)

        #  Occupation temperature
        label = myQtW.QLabel("Occupation mixing")
        self._occmix = myQtW.QLineEdit()
        validator = QtGui.QDoubleValidator(self._occmix)
        self._occmix.setValidator(validator)
        self._occmix.setText("1.0")
        form_layout.addRow(label, self._occmix)

        # Machines
        group_box = myQtW.QGroupBox("Computer and project name")
        self._layout.addWidget(group_box)

        form_layout = myQtW.QFormLayout()
        group_box.setLayout(form_layout)

        label = myQtW.QLabel("Machine")
        self._machine = myQtW.QComboBox()
        self._machine.addItems(["Summit", "bluewater"])
        form_layout.addRow(label, self._machine)
        label = myQtW.QLabel("Project name")

        self._projname = myQtW.QLineEdit()
        validator = QtGui.QDoubleValidator(self._projname)
        self._projname.setValidator(validator)
        self._projname.setText("CHP107")
        form_layout.addRow(label, self._projname)

        label = myQtW.QLabel("queue type")
        self._queue = myQtW.QComboBox()
        self._queue.addItems(["batch", "debug"])
        form_layout.addRow(label, self._queue)

    def state(self):
        """
        @return A dictionary containing the widget state.
        """
        input_setup_lines = (
            "# description of run\n"
            'description = "' + self._description.text() + '"\n\n'
            "# Uncommenting this will print out more information\n"
            '#verbose="true"\n\n'
            "# In most cases LCAO or Restart but certain special scenarios\n"
            "# may require a Random or Modified LCAO start\n"
            '#start_mode="LCAO Start"\n'
            '#start_mode="Random Start"\n'
            '#start_mode="Modified LCAO Start"\n'
            '#start_mode="Restart From File"\n'
            'start_mode = "' + self._start_mode.currentText() + '"\n\n'
            "# This is not an exhaustive list of options but does\n"
            "# contain the most frequently used ones.\n"
            '#calculation_mode="Quench Electrons"\n'
            '#calculation_mode="Relax Structure"\n'
            '#calculation_mode="Constant Volume And Energy"\n'
            '#calculation_mode="Constant Temperature And Energy"\n'
            '#calculation_mode="Band Structure Only"\n'
            'calculation_mode = "' + self._calculation_mode.currentText() + '"\n'
            'relax_method ="' + self._relaxmethod.currentText() + '"\n\n'
            "# Most pseudopotentials specify the exchange correlation type they\n"
            "# were generated with and the default value of AUTO_XC means that\n"
            "# the type specified in the pseudopotential is what RMG will use. That\n"
            "# can be overridden by specifying a value here. For a full list of\n"
            "# the available types look in the source distribution at the file\n"
            "# Headers/InputOpts.h around line 146.\n"
            'exchange_correlation_type = "' + self._ecf.currentText() + '"\n\n'
            "# RMG supports the following lattice types (Hexagonal at gamma-only)\n"
            '#bravais_lattice_type="Cubic Primitive"\n'
            '#bravais_lattice_type="Orthorhombic Primitive"\n'
            '#bravais_lattice_type="Hexagonal Primitive"\n'
            'bravais_lattice_type = "' + self._brav.currentText() + '"\n'
            'system_charge = "' + str(self._charge.text()) + '"\n'
        )

        input_occupation_lines = (
            "# RMG supports several different ways of specifying orbital occupations.\n"
            "# For a spin polarized system one may specify the occupations for up and\n"
            "# down separately. In the case of a non-zero electronic temperature these\n"
            "# will be adjusted as the calculation proceeds. For a non-spin polarized\n"
            "# calculation look at one of the other examples.\n"
            '#occupations_type = "Fixed"\n'
            '#occupations_type = "Fermi Dirac"\n'
            '#occupations_type = "MethfesselPaxton"\n'
            'occupations_type = "' + self._occ.currentText() + '"\n'
            'occupation_electron_temperature_eV = "' + str(self._occtem.text()) + '"\n'
            'occupation_number_mixing = "' + str(self._occmix.text()) + '"\n\n'
        )
        input_units_lines = (
            "#  length unit in Bohr or Angstrom\n"
            "# Default is Bohr\n"
            '#crds_units="Bohr"\n'
            '#crds_units="Angstrom"\n'
            'crds_units = "' + self._lengthunit.currentText() + '"\n'
            '#atomic_coordinate_type="Cell Relative"\n'
            '#atomic_coordinate_type="Absolute"\n'
            'atomic_coordinate_type = "' + self._atomcoor.currentText() + '"\n\n'
        )

        state = {
            "input_setup_lines": input_setup_lines,
            "input_units_lines": input_units_lines,
            "input_occupation_lines": input_occupation_lines,
            "length_units": str(self._lengthunit.currentText()) + "\n",
        }
        return state

    # end of state(self):
