# written by Wenchang Lu at NCSU
#
# Update Jun 24, 2023:
# - Update by Jackson Burns for Python 3

from PyQt5 import QtWidgets as myQtW

_NCSMURG_ADDON_PATH = "/AddOns/NCSURMG/"


class IOcontrol(myQtW.QWidget):
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
        group_box = myQtW.QGroupBox("IOcontrol")
        self._layout.addWidget(group_box)

        form_layout = myQtW.QFormLayout()
        group_box.setLayout(form_layout)

        # controle eigen value write

        label = myQtW.QLabel("Write eigvals period ")
        self._eigp = myQtW.QSpinBox()
        self._eigp.setValue(10)
        form_layout.addRow(label, self._eigp)

        # controle wavefunction  write

        label = myQtW.QLabel("wave function output after this md steps  ")
        self._wavemd = myQtW.QSpinBox()
        self._wavemd.setValue(10)
        form_layout.addRow(label, self._wavemd)

        # input wave files
        label = myQtW.QLabel("Input wave function file")
        self._inputwave = myQtW.QLineEdit()
        form_layout.addRow(label, self._inputwave)
        self._inputwave.setText("Wave/wave")

        # output wave files
        label = myQtW.QLabel("output wave function file")
        self._outputwave = myQtW.QLineEdit()
        form_layout.addRow(label, self._outputwave)
        self._outputwave.setText("Wave/wave")

    def state(self):
        """
        @return A dictionary containing the widget state.
        """
        input_io_lines = (
            "\n#  **** IO control  ****\n\n"
            'write_eigvals_period = "' + str(self._eigp.text()) + '"\n'
            'md_steps_til_write_waves = "' + str(self._wavemd.text()) + '"\n'
            'input_wave_function_file = "' + self._inputwave.text() + '"\n'
            'output_wave_function_file = "' + self._outputwave.text() + '"\n'
        )
        state = {"input_io_lines": input_io_lines}
        return state

    # end of state(self):
