# PP selector for NCSURMG plugin
# originally written by S.V. Barabash from QuantumWise
# modified by Wenchang Lu at NCSU
#
# Update Jun 24, 2023:
# - Update by Jackson Burns for Python 3

# defaults: ideally should be customizeable
# but keeping them fixed as requested by Wenchang

import warnings

from PyQt5 import QtCore, QtGui
from PyQt5 import QtWidgets
from PyQt5 import QtWidgets as myQtW

PP_PREFIX = "../"
PP_SUFFIX = ".pp"


class species(QtWidgets.QWidget):
    """
    Widget for selecting pseudopotentials.
    """

    pp_prefix = PP_PREFIX
    pp_suffix = PP_SUFFIX

    def __init__(self, parent=None):
        """
        Constructor.
        @param parent : The parent widget.
        """

        try:
            # Initialize base class
            myQtW.QWidget.__init__(self, parent)

            # list of present elements - to ensure the correct order
            # do not just use keys of the dict_pp_dirs
            self.elements = []

            # maximum number of allowed elements:
            self.max_n = 20

            #  __init__ : set up GUI:

            # Main layout
            self._layout = myQtW.QVBoxLayout()
            self.setLayout(self._layout)

            # Setup Groupbox
            group_box = myQtW.QGroupBox(
                "    Pseudopotential files       # of Orbitals      orbital radius  Hubbard_U"
            )
            self._layout.addWidget(group_box)

            form_layout = myQtW.QFormLayout()
            group_box.setLayout(form_layout)

            self.error_label = myQtW.QLabel("(No valid configuration defined)")
            form_layout.addRow(self.error_label)

            # create an array of max_n copies of PP selectors
            self.pp_layouts = list(range(self.max_n))
            self.pp_labels = list(range(self.max_n))
            self.pp_labels1 = list(range(self.max_n))
            self.pp_labels2 = list(range(self.max_n))
            self.pp_lines = list(range(self.max_n))
            self.pp_buttons = list(range(self.max_n))
            self.num_orbital = list(range(self.max_n))
            self.orbital_radius = list(range(self.max_n))
            self.hubbard_u = list(range(self.max_n))

            # first define slot-like versions of the
            # selectPPDir() function, to be used for the
            # connections below
            def createPPConnect(i_pp):
                """
                internal lambda-creating function to connect
                parameter-free signal from radio buttons for PP variants to
                a parameter-accepting class method
                """
                return lambda: self.selectPPDir(i_pp)

            for i_pp in range(self.max_n):
                self.pp_layouts[i_pp] = myQtW.QHBoxLayout()

                self.pp_labels[i_pp] = myQtW.QLabel("undefined")
                self.pp_labels1[i_pp] = myQtW.QLabel(" ")
                self.pp_labels2[i_pp] = myQtW.QLabel(" ")
                self.pp_lines[i_pp] = myQtW.QLineEdit("undefined")
                self.pp_buttons[i_pp] = myQtW.QPushButton("Browse...")
                self.num_orbital[i_pp] = myQtW.QSpinBox()
                self.num_orbital[i_pp].setValue(6)

                self.orbital_radius[i_pp] = myQtW.QLineEdit()
                validator = QtGui.QDoubleValidator(self.orbital_radius[i_pp])
                self.orbital_radius[i_pp].setValidator(validator)
                self.orbital_radius[i_pp].setText("8.0")

                self.hubbard_u[i_pp] = myQtW.QLineEdit()
                validator = QtGui.QDoubleValidator(self.hubbard_u[i_pp])
                self.hubbard_u[i_pp].setValidator(validator)
                self.hubbard_u[i_pp].setText("0.0")

                self.pp_layouts[i_pp].addWidget(self.pp_labels[i_pp])
                self.pp_layouts[i_pp].addWidget(self.pp_lines[i_pp])
                self.pp_layouts[i_pp].addWidget(self.pp_buttons[i_pp])
                self.pp_layouts[i_pp].addWidget(self.pp_labels1[i_pp])
                self.pp_layouts[i_pp].addWidget(self.num_orbital[i_pp])
                self.pp_layouts[i_pp].addWidget(self.pp_labels2[i_pp])
                self.pp_layouts[i_pp].addWidget(self.orbital_radius[i_pp])
                self.pp_layouts[i_pp].addWidget(self.hubbard_u[i_pp])

                self.pp_lines[i_pp].setToolTip("edit the pseudopotential path")
                self.pp_buttons[i_pp].setToolTip("select an existing directory")

                self.pp_buttons[i_pp].clicked.connect(createPPConnect(i_pp))

                form_layout.addRow(self.pp_layouts[i_pp])

            self.setElements(None)

        except Exception as e:
            print("error in species", e)

    # end of  __init__

    def setElements(self, configuration):
        """
        Method called externally (by other modules)
        when the set of the available
        elements changes (i.e. when a configuration changes).
        """
        if configuration is None:
            self.elements = []
            self.error_label.setVisible(True)
        else:
            # valid configuration is supplied, build the list of element types:
            _elements = configuration.conf.elements
            _element_types = list(set(_elements))

            if len(_element_types) > 0 and len(_element_types) < self.max_n:
                self.elements = _element_types
                self.error_label.setVisible(False)

            else:
                self.elements = []
                warnings.warn(
                    "The configuration may not have more than "
                    + str(self.max_n)
                    + " elements or no elements at all!"
                )

        # make the necessary pp controls visible and properly initialized:
        for i_pp in range(len(self.elements)):
            self.pp_labels[i_pp].setText(self.elements[i_pp])
            _default_name = ""
            self.pp_lines[i_pp].setText(_default_name)
            self.pp_labels[i_pp].setVisible(True)
            self.pp_labels1[i_pp].setVisible(True)
            self.pp_labels2[i_pp].setVisible(True)
            self.pp_lines[i_pp].setVisible(True)
            self.pp_buttons[i_pp].setVisible(True)
            self.num_orbital[i_pp].setVisible(True)
            self.orbital_radius[i_pp].setVisible(True)
            self.hubbard_u[i_pp].setVisible(True)

        # make the other controls invisible:
        for i_pp in range(len(self.elements), self.max_n):
            self.pp_labels[i_pp].setVisible(False)
            self.pp_labels1[i_pp].setVisible(False)
            self.pp_labels2[i_pp].setVisible(False)
            self.pp_lines[i_pp].setVisible(False)
            self.pp_buttons[i_pp].setVisible(False)
            self.num_orbital[i_pp].setVisible(False)
            self.orbital_radius[i_pp].setVisible(False)
            self.hubbard_u[i_pp].setVisible(False)

    def selectPPDir(self, i_pp):
        """
        The method that is called when the user selects to
        browse to a directory for the i_pp'th pseudopotential
        """
        try:
            dialog = myQtW.QFileDialog(self)
            #            dialog.setFileMode(myQtW.QFileDialog.Directory)
            if QtCore.QDir(self.pp_lines[i_pp].text()).exists():
                dialog.selectFile(self.pp_lines[i_pp].text())
            if dialog.exec_():
                filename = dialog.selectedFiles()[0]
                self.pp_lines[i_pp].setText(filename)
        except Exception as e:
            print(e)
            warnings.warn(
                "Error selecting pseudopotential id = {:d} ".format(
                    i_pp,
                )
            )

    def state(self):
        """
        @return A dictionary containing the widget state.
        """

        species_pp = ""
        for i in range(len(self.elements)):
            if str(self.pp_lines[i].text()) != "":
                species_pp += (
                    " "
                    + self.elements[i]
                    + "    "
                    + str(self.pp_lines[i].text())
                    + "\n"
                )

        species_lines = ""
        if species_pp == "":
            species_lines = "\n# Using internal Pseudopotential **** \n\n"
        else:
            species_lines += 'pseudopotential="\n'
            species_lines += species_pp
            species_lines += '"\n'

        hubbard_line = ""
        for i in range(len(self.elements)):
            print(self.hubbard_u[i].text())
            if self.hubbard_u[i].text() != "0.0":
                hubbard_line += (
                    self.elements[i] + "   " + str(self.hubbard_u[i].text()) + "\n"
                )
                ldaU_radius = self.orbital_radius[i].text()
        if hubbard_line != "":
            species_lines += 'ldaU_mode = "Simple"' + "\n"
            species_lines += 'ldaU_radius = "' + ldaU_radius + '"\n'
            species_lines += 'Hubbard_U="\n'
            species_lines += hubbard_line
            species_lines += '"\n'

        state = {"input_species_lines": species_lines}

        return state

    # end of state(self):
