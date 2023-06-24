# written by Wenchang Lu at NCSU
#
# Update Jun 24, 2023:
# - Update by Jackson Burns for Python 3


from PyQt5 import QtGui
from PyQt5 import QtWidgets as myQtW

_NCSMURG_ADDON_PATH = "/AddOns/NCSURMG/"


class Mdscf(myQtW.QWidget):
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

        # __init__ : set up GUI:

        # Main layout
        self._layout = myQtW.QVBoxLayout()
        self.setLayout(self._layout)

        # Setup Groupbox
        group_box = myQtW.QGroupBox("SCF")
        self._layout.addWidget(group_box)

        form_layout = myQtW.QFormLayout()
        group_box.setLayout(form_layout)

        # max SCF steps

        label = myQtW.QLabel("Max scf steps for RMG ")
        self._maxscf = myQtW.QLineEdit()
        self._maxscf.setText("50")
        form_layout.addRow(label, self._maxscf)

        # RMS convergence critierion
        label = myQtW.QLabel("RMS Convergence Criterion")
        self._rms = myQtW.QLineEdit()
        validator = QtGui.QDoubleValidator(self._rms)
        self._rms.setValidator(validator)
        form_layout.addRow(label, self._rms)
        self._rms.setText("1e-7")

        label = myQtW.QLabel("energy Convergence Criterion")
        self._rms_energy = myQtW.QLineEdit()
        validator = QtGui.QDoubleValidator(self._rms)
        self._rms_energy.setValidator(validator)
        form_layout.addRow(label, self._rms_energy)
        self._rms_energy.setText("1e-9")

        # Setup Groupbox
        group_box = myQtW.QGroupBox("Relax/MD")
        self._layout.addWidget(group_box)

        form_layout = myQtW.QFormLayout()
        group_box.setLayout(form_layout)
        # max MD steps

        label = myQtW.QLabel("Max MD/FastRelax steps ")
        self._maxmd = myQtW.QLineEdit()
        self._maxmd.setText("10")
        form_layout.addRow(label, self._maxmd)

        # force convergence critierion
        label = myQtW.QLabel("Max force (Ha/au) < ")
        self._forcemax = myQtW.QLineEdit()
        validator = QtGui.QDoubleValidator(self._forcemax)
        self._forcemax.setValidator(validator)
        form_layout.addRow(label, self._forcemax)
        self._forcemax.setText("1e-3")

        # ionic time step
        label = myQtW.QLabel("ionic time step")
        self._ionstep = myQtW.QLineEdit()
        validator = QtGui.QDoubleValidator(self._ionstep)
        self._ionstep.setValidator(validator)
        form_layout.addRow(label, self._ionstep)
        self._ionstep.setText("45")

        # relax_dynamic_timestep
        label = myQtW.QLabel("time step changed dynamicall")
        self._rdt = myQtW.QComboBox()
        self._rdt.addItems(["true", "false"])
        form_layout.addRow(label, self._rdt)

        # Setup Groupbox
        group_box = myQtW.QGroupBox("Mixing")
        self._layout.addWidget(group_box)
        form_layout = myQtW.QFormLayout()
        group_box.setLayout(form_layout)

        # charge density mixing parameter
        label = myQtW.QLabel("Charge Density mixing")
        self._qmix = myQtW.QLineEdit()
        validator = QtGui.QDoubleValidator(self._qmix)
        self._qmix.setValidator(validator)
        form_layout.addRow(label, self._qmix)
        self._qmix.setText("0.7")

        # projector mixing parameter
        label = myQtW.QLabel("Projector mixing")
        self._pmix = myQtW.QLineEdit()
        validator = QtGui.QDoubleValidator(self._pmix)
        self._pmix.setValidator(validator)
        # form_layout.addRow(label, self._pmix)
        self._pmix.setText("0.1")

        #  charge density mixing method
        label = myQtW.QLabel("charge density mixing method")
        self._mixmethod = myQtW.QComboBox()
        self._mixmethod.addItems(["Broyden", "Pulay", "Linear"])
        form_layout.addRow(label, self._mixmethod)

        Hlayout = myQtW.QHBoxLayout()
        # Pulay mixing histroy
        label = myQtW.QLabel("  Pulay Order ")
        self._pulayorder = myQtW.QSpinBox()
        self._pulayorder.setValue(5)
        Hlayout.addWidget(label)
        Hlayout.addWidget(self._pulayorder)

        # Pulay mixing scale
        label = myQtW.QLabel(" scale (beta)  ")
        self._pulaybeta = myQtW.QLineEdit()
        self._pulaybeta.setText("0.5")
        Hlayout.addWidget(label)
        Hlayout.addWidget(self._pulaybeta)

        # Pulay mixing refresh
        label = myQtW.QLabel(" refresh steps  ")
        self._pulayrefresh = myQtW.QLineEdit()
        self._pulayrefresh.setText("100")
        Hlayout.addWidget(label)
        Hlayout.addWidget(self._pulayrefresh)
        form_layout.addRow(Hlayout)

    def state(self):
        """
        @return A dictionary containing the widget state.
        """

        input_mdscf_lines = (
            "\n# **** MD/Relax controls  ****\n  \n"
            'max_md_steps = "' + str(self._maxmd.text()) + '"\n'
            'relax_max_force = "' + str(self._forcemax.text()) + '"\n'
            'ionic_time_step = "' + str(self._ionstep.text()) + '"\n'
            'relax_dynamic_timestep = "{:s}"\n'
            "\n# **** Mixing controls **** \n  \n"
            "# RMG supports Broyden, Pulay and Linear mixing\n"
            "# When the davidson Kohn-Sham solver is selected Broyden or\n"
            "# Pulay are preferred. For the multigrid solver Linear with\n"
            "# potential acceleration is often (but not always) the best\n"
            "# choice.\n"
            '#charge_mixing_type = "Broyden"\n'
            '#charge_mixing_type = "Pulay"\n'
            '#charge_mixing_type = "Linear"\n'
            'charge_mixing_type = "{:s}"\n'
            'charge_density_mixing = "' + str(self._qmix.text()) + '"\n'
            'charge_pulay_order = "{:s}"\n'
            'charge_pulay_scale = "{:s}"\n'
            'charge_pulay_refresh = "{:s}"\n'
            "\n# **** SCF controls ****  \n  \n"
            'max_scf_steps = "' + str(self._maxscf.text()) + '"\n'
            'rms_convergence_criterion = "' + str(self._rms.text()) + '"\n'
            'energy_convergence_criterion = "' + str(self._rms_energy.text()) + '"\n'
        ).format(
            self._rdt.currentText(),
            self._mixmethod.currentText(),
            self._pulayorder.text(),
            self._pulaybeta.text(),
            self._pulayrefresh.text(),
        )
        state = {"input_mdscf_lines": input_mdscf_lines}
        return state

    # end of state(self):
