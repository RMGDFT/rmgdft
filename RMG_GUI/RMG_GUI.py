# Wenchang Lu at NCSU
#
# Update Jun 24, 2023:
# - Update by Jackson Burns for Python 3

import datetime
import os
import sys
import warnings

from PyQt5 import QtWidgets as myQtW

from Configuration import Configuration
from default_input_para import default_input_para
from Grids import Grids
from IOcontrol import IOcontrol
from Mdscf import Mdscf
from Misc import Misc
from Setup import Setup
from species import species
from write_out_jobfiles import write_out_jobfiles


class RMG_GUI(myQtW.QTabWidget):
    """
    Class used for representing a NCSURMG scripter.
    """

    def __init__(self):
        super(RMG_GUI, self).__init__()

        self.initUI()

    def initUI(self):
        # Add the widgets
        self._setup = Setup()
        self._species = species()
        self._misc = Misc()
        self._grids = Grids()
        self._mdscf = Mdscf()
        self._io = IOcontrol()
        self._configuration = Configuration()
        self._default_input = default_input_para()
        self._widgets = [
            self._setup,
            self._misc,
            self._mdscf,
            self._io,
            self._species,
            self._grids,
            self._default_input,
        ]
        # Main layout
        layout = myQtW.QGridLayout()
        layout.setSpacing(10)
        self.setLayout(layout)

        # Setup Groupbox
        savebtn = myQtW.QPushButton("Save")
        self.savedir = myQtW.QLineEdit(os.getcwd())
        choosedir = myQtW.QPushButton("...")
        choosedir.clicked.connect(self.selectdir)
        savebtn.clicked.connect(self.save)
        layout.addWidget(savebtn, 1, 0, 1, 1)
        layout.addWidget(self.savedir, 1, 1, 1, 4)
        layout.addWidget(choosedir, 1, 5, 1, 1)
        form_layout = myQtW.QTabWidget()
        layout.addWidget(form_layout, 2, 0, 1, 10)
        form_layout.addTab(self._setup, self.tr("Setup"))
        form_layout.addTab(self._mdscf, self.tr("MD SCF"))
        form_layout.addTab(self._grids, self.tr("Grids"))
        form_layout.addTab(self._species, self.tr("Species"))
        form_layout.addTab(self._misc, self.tr("Misc"))
        form_layout.addTab(self._io, self.tr("IO"))
        form_layout.addTab(self._configuration, self.tr("Configuration"))
        self._species.setElements(self._configuration)
        self._grids.lattparameters(self._configuration, self._misc)
        self._configuration.CONF_CHANGED.connect(self.configurationChanged)
        self.setGeometry(2000, 2000, 750, 750)
        self.setWindowTitle("NCSU RMG GUI")
        self.show()

    # end of __init__
    def update(self):
        if self._configuration.input_para != {}:
            calc_modes = [
                self._setup._calculation_mode.itemText(i)
                for i in range(self._setup._calculation_mode.count())
            ]
        calc_mode = self._configuration.input_para["calculation_mode"]
        ip = self._configuration.input_para
        self._setup._calculation_mode.setCurrentIndex(
            calc_modes.index(calc_mode),
        )
        self._setup._description.setText(ip["description"])
        rlx_methods = [
            self._setup._relaxmethod.itemText(i)
            for i in range(self._setup._relaxmethod.count())
        ]
        rlx_method = ip["relax_method"]
        self._setup._relaxmethod.setCurrentIndex(rlx_methods.index(rlx_method))
        ecfs = [self._setup._ecf.itemText(i) for i in range(self._setup._ecf.count())]
        self._setup._ecf.setCurrentIndex(
            ecfs.index(ip["exchange_correlation_type"]),
        )
        bravs = [
            self._setup._brav.itemText(i) for i in range(self._setup._brav.count())
        ]
        self._setup._brav.setCurrentIndex(bravs.index(ip["bravais_lattice_type"]))
        self._setup._charge.setText(ip["system_charge"])
        self._setup._occtem.setText(ip["occupation_electron_temperature_eV"])
        self._setup._occmix.setText(ip["occupation_number_mixing"])
        self._setup._lengthunit.setText(ip["length_units"])
        self._setup._atomcoor.setText(ip["atomic_coordinate_type"])
        self._mdscf._maxmd.setText(ip["max_md_steps"])
        self._mdscf._forcemax.setText(ip["relax_max_force"])
        self._mdscf._ionstep.setText(ip["ionic_time_step"])
        self._mdscf._rdt.setCurrentIndex(
            self._mdscf._rdt.findText(ip["relax_dynamix_timestep"])
        )
        self._mdscf._qmix.setText(ip["charge_desity_mixing"])
        self._mdscf._pmix.setText(ip["projector_mixing"])
        self._mdscf._mixmethod.setCurrentIndex(
            self._mdscf._mixmethod.findText(ip["charge_mixing_type"])
        )
        self._mdscf._pulayorder.setText(ip["charge_pulay_order"])
        self._mdscf._pulaybeta.setText(ip["charge_pulay_scale"])
        self._mdscf._pulayrefresh.setText(ip["charge_pulay_refresh"])
        self._mdscf._maxscf.setText(ip["max_scf_steps"])
        self._mdscf._rms.setText(ip["rms_convergence_criterion"])
        self._misc._khlevel.setText(ip["kohn_sham_mg_levels"])
        self._misc._poissonlevel.setText(ip["poisson_mg_levels"])
        self._misc._hkstep.setText(ip["kohn_sham_time_step"])
        self._misc._poissonstep.setText(ip["poisson_time_step"])
        self._io._eigp.setText(ip["write_eigvals_period"])
        self._io._wavemd.setText(ip["md_steps_till_write_waves"])
        self._io._inputwave.setText(ip["input_wave_function_file"])
        self._io._outputwave.setText(ip["output_wave_function_file"])

    def configurationChanged(self):
        """
        Called automatically when a configuration is dropped on the tool.
        @param configuration : The new configuration.
        """
        # Set the configuration

        # Move focus to the configuration tab
        self.setCurrentWidget(self._configuration)
        self.setCurrentWidget(self._grids)
        self._species.setElements(self._configuration)
        self._grids.lattparameters(self._configuration, self._misc)

    def save(self):
        """
        @Make the files
        """
        directory = self.savedir.text()
        os.chdir(directory)
        self.setCurrentWidget(self._setup)
        self.setCurrentWidget(self._grids)
        self.setCurrentWidget(self._species)
        self.setCurrentWidget(self._configuration)
        _mystate = self.state()
        configuration = self._configuration
        _elements = []
        _positions_line = "# undefined configuration\n"
        _elements = configuration.conf.elements

        # write positions:
        _lattice_line = (
            "#  **** Lattice constants **** \n\n"
            'a_length =" ' + str(configuration.conf.lattice[0]) + '"\n'
            'b_length =" ' + str(configuration.conf.lattice[1]) + '"\n'
            'c_length =" ' + str(configuration.conf.lattice[2]) + '"\n'
        )

        _positions_line = ""
        zipped_coordinates = zip(_elements, configuration.conf.coords)
        _positions_line += 'atoms=\n"'
        for name, v in zipped_coordinates:
            _positions_line += (
                " {:s}     {:.12e}    {:.12e}    {:.12e}      1   0.0\n".format(
                    name, v[0], v[1], v[2]
                )
            )
        _positions_line += '"\n'
        _positions_line += _lattice_line

        common_lines = """
# RMGDFT Input File
# usage: /path/to/rmg rmg.input
#
# This input file was automatically generated by RMG_GUI.py
# at this timestamp: {:s}
#
# It may need to be edited further to be compatible with the
# current version of RMG-DFT.

""".format(
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )
        common_lines += _mystate["input_setup_lines"]
        common_lines += _mystate["input_grids_lines"]
        common_lines += _mystate["input_units_lines"]
        common_lines += _mystate["input_occupation_lines"]
        common_lines += _mystate["input_misc_lines"]
        common_lines += _mystate["input_io_lines"]
        common_lines += _mystate["input_species_lines"]
        common_lines += _mystate["input_mdscf_lines"]
        common_lines += _positions_line

        filename = "/rmg.input"
        if self._setup._calculation_mode.currentText() == "Band Structure Only":
            filename += ".band"

        print(
            "Attempting to write current configuration to",
            directory + filename,
        )
        with open(directory + filename, "w") as input_file:
            input_file.write(common_lines)
            print("Wrote RMG input file to", directory + filename)

        try:
            write_out_jobfiles(configuration, self._setup, self._grids)
        except Exception as e:
            warnings.warn("Unable to write job to slurm file,", str(e))

    def state(self):
        """
        @return The state of all widgets in a state dictionary.
        """
        state = {}
        for widget in self._widgets:
            state.update(widget.state())
        return state

    def title(self):
        """
        @return The title of the plugin.
        """
        return "QW-RMG NCSU"

    def selectdir(self):
        directory = myQtW.QFileDialog.getExistingDirectory(self)
        self.savedir.setText(directory)


def main():
    app = myQtW.QApplication(sys.argv)
    ex = RMG_GUI()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
