# Wenchang Lu at NCSU

from PyQt4 import QtCore, QtGui

import os
import sys

from Setup import Setup
#from RMG_para import rmg_para
from lattice import lattice
#from write_out_LCR import write_out_LCR
#from write_out_others import write_out_others
from write_out_jobfiles import write_out_jobfiles
#from position_orbital import position_orbital
#from position_orbital_3part import position_orbital_3part
#from position_orbital_negf import position_orbital_negf
#from position_orbital_center_ON import position_orbital_center_ON
#from position_orbital_center_ON_nop import position_orbital_center_ON_nop
from default_input_para import default_input_para
from Misc import Misc
from Grids import Grids
from IOcontrol import IOcontrol
from species import species
from Mdscf import Mdscf
from Configuration import Configuration

from distutils.sysconfig import get_python_lib


    
########################### begin class definition



class RMG_GUI(QtGui.QTabWidget):
    """
       Class used for representing a NCSURMG scripter.
    """
       
    try:
        
        def __init__(self):
            super(RMG_GUI, self).__init__()
        
            self.initUI()
        
        def initUI(self):               
        

            # Add the widgets
            self._setup = Setup()
#            self._rmg_para = rmg_para()
            self._species=species()
            self._misc = Misc()
            self._grids = Grids()
            self._mdscf = Mdscf()
            self._io = IOcontrol()
            self._configuration = Configuration()
            self._default_input = default_input_para()

            #self._widgets = [self._setup, self._misc, self._mdscf, self._io, 
            #                 self._species, self._grids, self._default_input,
            #                 self._configuration, self._negf_para]
            self._widgets = [self._setup, self._misc, self._mdscf, self._io,
                             self._species, self._grids, self._default_input]
			    #,self._rmg_para]#,self._configuration]
    # Main layout
            #layout = QtGui.QVBoxLayout()
            layout = QtGui.QGridLayout()
            layout.setSpacing(10)
            self.setLayout(layout)

        # Setup Groupbox
            savebtn = QtGui.QPushButton('Save')
            self.savedir = QtGui.QLineEdit(os.getcwd())
            choosedir = QtGui.QPushButton('...')
            choosedir.clicked.connect(self.selectdir)
            savebtn.clicked.connect(self.save)

            layout.addWidget(savebtn, 1, 0, 1, 1)
            layout.addWidget(self.savedir, 1, 1, 1, 4)
            layout.addWidget(choosedir, 1, 5, 1, 1)

            form_layout = QtGui.QTabWidget()
            layout.addWidget(form_layout, 2, 0, 1, 10)

            #form_layout = QtGui.QTabWidget(group_box)


            
            form_layout.addTab(self._setup,         self.tr('Setup'))
            form_layout.addTab(self._mdscf,         self.tr('MD SCF'))
            form_layout.addTab(self._grids,          self.tr('Grids'))
#            form_layout.addTab(self._rmg_para,          self.tr('RMGctrl'))
            form_layout.addTab(self._species,         self.tr('Species'))
            form_layout.addTab(self._misc,          self.tr('Misc'))
            form_layout.addTab(self._io,            self.tr('IO'))
            form_layout.addTab(self._configuration, self.tr('Configuration'))

            
            self._species.setElements(self._configuration)
            self._grids.lattparameters(self._configuration, self._misc)
            #form_layout.currentChanged.connect(self.configurationChanged)
#            self.connect(form_layout, QtCore.SIGNAL("currenChanged(int)"), self.configurationChanged)
            self._configuration.CONF_CHANGED.connect(self.configurationChanged)
#     self._configuration.button.clicked.connect(self.update)
            
            self.setGeometry(2000, 2000, 750, 750)
            self.setWindowTitle('NCSU RMG GUI')    
            self.show()

        ######## end of __init__
	def update(self):
	    if self._configuration.input_para!={}:
		#print 'up'
	        #print self._configuration.input_para

	    	calc_modes=[self._setup._calculation_mode.itemText(i) for i in range(self._setup._calculation_mode.count())]
		calc_mode=self._configuration.input_para['calculation_mode']
		ip=self._configuration.input_para
		self._setup._calculation_mode.setCurrentIndex(calc_modes.index(calc_mode))
		self._setup._description.setText(ip['description'])
	 	rlx_methods=[self._setup._relaxmethod.itemText(i) for i in range(self._setup._relaxmethod.count())]
	   	rlx_method=ip['relax_method']
		self._setup._relaxmethod.setCurrentIndex(rlx_methods.index(rlx_method))
		ecfs=[self._setup._ecf.itemText(i) for i in range(self._setup._ecf.count())]
		self._setup._ecf.setCurrentIndex(ecfs.index(ip['exchange_correlation_type']))
		bravs=[self._setup._brav.itemText(i) for i in range(self._setup._brav.count())]
		self._setup._brav.setCurrentIndex(bravs.index(ip['bravais_lattice_type']))


		self._setup._charge.setText(ip['system_charge'])

#		iccs=[self._setup._occ.itemText(i) for i in range(self._setup._occ.count())]
#self._setup._occ.setCurrentIndex(self._setup._occ.findText(ip['occupation_type']))
		self._setup._occtem.setText(ip['occupation_electron_temperature_eV'])
		self._setup._occmix.setText(ip['occupation_number_mixing'])
		self._setup._lengthunit.setText(ip['length_units'])
		self._setup._atomcoor.setText(ip['atomic_coordinate_type'])



		self._mdscf._maxmd.setText(ip['max_md_steps'])
		self._mdscf._forcemax.setText(ip['relax_max_force'])
		self._mdscf._ionstep.setText(ip['ionic_time_step'])
#		rdts=[self._mdscf._rdt.itemText(i) for i in range(self._mdscf._rdt.count())]
#		self._mdscf._rdt.setCurrentIndex
		self._mdscf._rdt.setCurrentIndex(self._mdscf._rdt.findText(ip['relax_dynamix_timestep']))
		self._mdscf._qmix.setText(ip['charge_desity_mixing'])
		self._mdscf._pmix.setText(ip['projector_mixing'])
		self._mdscf._mixmethod.setCurrentIndex(self._mdscf._mixmethod.findText(ip['charge_mixing_type']))
		self._mdscf._pulayorder.setText(ip['charge_pulay_order'])
		self._mdscf._pulaybeta.setText(ip['charge_pulay_scale'])
		self._mdscf._pulayrefresh.setText(ip['charge_pulay_refresh'])		
		self._mdscf._maxscf.setText(ip['max_scf_steps'])
		self._mdscf._rms.setText(ip['rms_convergence_criterion'])



		self._misc._khlevel.setText(ip['kohn_sham_mg_levels'])
		self._misc._poissonlevel.setText(ip['poisson_mg_levels'])
		self._misc._hkstep.setText(ip['kohn_sham_time_step'])
		self._misc._poissonstep.setText(ip['poisson_time_step'])

		

		self._io._eigp.setText(ip['write_eigvals_period'])
		self._io._wavemd.setText(ip['md_steps_till_write_waves'])
		self._io._inputwave.setText(ip['input_wave_function_file'])
		self._io._outputwave.setText(ip['output_wave_function_file'])
		




        def configurationChanged(self):
            """
               Called automatically when a configuration is dropped on the tool.            
               @param configuration : The new configuration.
            """
            # Set the configuration

            # Move focus to the configuration tab
            self.setCurrentWidget(self._configuration)
#            self.setCurrentWidget(self._setup)
            self.setCurrentWidget(self._grids)
            self._species.setElements(self._configuration)
            self._grids.lattparameters(self._configuration, self._misc)
		#self.update()

	#	calc_modes=[
        #                    "Quench Electrons",
        #                    "Relax Structure",
        #                    "NEB Relax",
        #                    "Constant Volume And Energy",
        #                    "Constant Temperature And Energy",
        #                    "Constant Pressure And Energy(not implemented)",
        #                    "Constrained Fast Relax",
        #                    "Band Structure Only"
        #                     ]
#            self._grids.lattparameters(self._configuration)
#            self._grids.changeothergrid()
#            self._grids.get_nx_negf()
            
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
            #self.configurationChanged(self._configuration)
            #self.configurationChanged()

#            try:
            _mystate = self.state()
            
            configuration = self._configuration
            _format_vector = "  {0: 10.8f} {1: 10.8f} {2: 10.8f}\n"
            _format_with_comment = "  {0: 10.8f} {1: 10.8f} {2: 10.8f} {3}\n"
            _ncsu_atom_format = " %s     %.12e    %.12e    %.12e      1   0.0\n"


            _elements = []
            _element_types = []
            _element_count = []
            _positions_line="# undefined configuration\n"


            zipped_counts = {}
            # Test whether we have a configuration defined:


            #e generate types of atoms encountered:

            _elements = configuration.conf.elements
            _element_types = list(set(_elements))
            _element_count = [ _elements.count(ii) for ii in _element_types]

            
            # write positions:

            _lattice_line = "#  **** Lattice constants **** \n\n"
	    _lattice_line +='a_length =" '+str(configuration.conf.lattice[0])+'"\n'
	    _lattice_line +='b_length =" '+str(configuration.conf.lattice[1])+'"\n'
	    _lattice_line +='c_length =" '+str(configuration.conf.lattice[2])+'"\n'

            
            _positions_line = ''
            zipped_coordinates = zip(_elements, configuration.conf.coords)
            _positions_line += 'atoms=\n"'
            for name,v in zipped_coordinates:
                _positions_line += ( _ncsu_atom_format % (name, v[0],v[1],v[2]) )
            _positions_line += '"\n'
            _positions_line +=_lattice_line


#define dicary for number of orbitals and radius
            num_orbital_dict = {}
            num_orbital = [self._species.num_orbital[i].value() for i in range(len(_element_types))]
            num_orbital_dict = dict(zip(_element_types, num_orbital))
            
            orbital_radius_dict = {}
            orbital_radius = [self._species.orbital_radius[i].text() for i in range(len(_element_types))]
            orbital_radius_dict = dict(zip(_element_types, orbital_radius))

            common_lines  =  _mystate['input_setup_lines']
            common_lines += _mystate['input_grids_lines']
            common_lines +=  _mystate['input_units_lines']
            common_lines +=  _mystate['input_occupation_lines']
            common_lines += _mystate['input_misc_lines']
            common_lines += _mystate['input_io_lines']
            common_lines += _mystate['input_species_lines']
#            common_lines += _mystate['input_species_lines']
            common_lines += _mystate['input_mdscf_lines']
#            common_lines += _mystate['input_conf_lines']
            common_lines += _positions_line

            filename = '/rmg.input'
            if(self._setup._calculation_mode.currentText() =='Band Structure Only'): filename = '/rmg.input.band'
            
            print filename 
            with open(directory + filename,'w') as input_file:
                input_file.write(common_lines)

#            except:
#               print "failed to save"
#            try:
#
##  Use the lattice parameter from left lead for the one atom calculation 
#                electrodes = configuration.conf[1]
#                latt = lattice(electrodes.lattice)
#
#
#                for i_pp in range(len(_element_types)):
#                    dir_name = _element_types[i_pp] + '-atom'
#                    oneatom_lines  = 'atoms=\n"\n'
#                    oneatom_lines += _element_types[i_pp] + '   0.0  0.0 0.0  1  1  '
#                    oneatom_lines += "%d" % num_orbital_dict[_element_types[i_pp]]
#                    oneatom_lines += '\n"\n\n'
#                    oneatom_lines += 'orbitals=\n"\n'
#                    oneatom_lines +=  str(self._species.num_orbital[i_pp].text())
#                    oneatom_lines +=  '  0.0  0.0  0.0 '
#                    oneatom_lines +=  orbital_radius_dict[_element_types[i_pp]]
#                    oneatom_lines +=  '  1  1\n"\n\n'
#                    oneatom_lines +=  'number_of_orbitals="'
#                    oneatom_lines +=  str(self._species.num_orbital[i_pp].text())
#                    oneatom_lines +=  '"\n\n'
#
#                
#                    if not os.path.exists(dir_name):
#                        os.mkdir(dir_name)
#                    with open(dir_name + '/input','w') as inc_file:
#                        inc_file.write(_mystate['input_oneatom_grid_lines'])
#                        inc_file.write(_mystate['input_mdscf_lines_ON'])
#                        inc_file.write(common_lines)
#                        inc_file.write(latt)
#                        inc_file.write(oneatom_lines)
#                        inc_file.write(_mystate['default_input_for_oneatom'])
#                        inc_file.write(_mystate['default_input_forON'])
##
##    write out order-n calculation input for lead1
#                        
#            except:
#               print "failed to save2"
#            try:
#                electrodes = configuration.conf[1]
#                latt = lattice(electrodes.lattice)
#                position = position_orbital(electrodes,
#                                   num_orbital_dict, orbital_radius_dict)
#                dir_name = 'lead1'
#                if not os.path.exists(dir_name): os.mkdir(dir_name)
#
#                with open(dir_name + '/input','w') as inc_file:
#                   inc_file.write('start_mode="FIREBALL Start"')
#                   inc_file.write(_mystate['input_grids_left_lines'])
#                   inc_file.write(_mystate['input_mdscf_lines_ON'])
#                   inc_file.write(common_lines)
#                   inc_file.write(latt)
#                   inc_file.write(position)
#                   inc_file.write(_mystate['default_input_forON'])
#
##
##    write out order-n calculation input for lead2
#                        
#            except:
#               print "failed to save2a"
#            try:
#                electrodes = configuration.conf[2]
#                latt = lattice(electrodes.lattice)
#                position = position_orbital(electrodes,
#                                   num_orbital_dict, orbital_radius_dict)
#                dir_name = 'lead2'
#                if not os.path.exists(dir_name): os.mkdir(dir_name)
#
#                with open(dir_name + '/input','w') as inc_file:
#                   inc_file.write('start_mode="FIREBALL Start"')
#                   inc_file.write(_mystate['input_grids_right_lines'])
#                   inc_file.write(_mystate['input_mdscf_lines_ON'])
#                   inc_file.write(common_lines)
#                   inc_file.write(latt)
#                   inc_file.write(position)
#                   inc_file.write(_mystate['default_input_forON'])
#
##
##    write out order-n calculation input for center part
#                        
#            except:
#               print "failed to save2b"
#            try:
#                center_ctrl = {}
#                if(self._grids.center_periodicity.isChecked()):
#                    center_ctrl = {'state_begin':0, 'ion_begin':0,'num_atoms':0}
#                    position = position_orbital_center_ON(configuration,
#                                   num_orbital_dict, orbital_radius_dict)
#                else:
#                    position,center_ctrl =position_orbital_center_ON_nop(configuration, 
#                                   num_orbital_dict, orbital_radius_dict)
#                dir_name = 'center'
#                if not os.path.exists(dir_name): os.mkdir(dir_name)
#
#                with open(dir_name + '/input','w') as inc_file:
#                   inc_file.write('start_mode="FIREBALL Start"')
#                   inc_file.write(_mystate['input_grids_center_lines'])
#                   inc_file.write(_mystate['input_mdscf_lines_ON'])
#                   inc_file.write(common_lines)
#                   inc_file.write(position)
#                   inc_file.write(_mystate['default_input_forON'])
#            except:
#               print "failed to save3"
#            try:
#
#                electrods = configuration.conf
#
#                conf1 = electrods[1]
#                conf2 = electrods[1]
#                conf3 = electrods[1]
#                position = position_orbital_3part(conf1,conf2,conf3,
#                                   num_orbital_dict, orbital_radius_dict)
#                
#                dir_name = '3lead_lead1'
#                if not os.path.exists(dir_name): os.mkdir(dir_name)
#                a1 = self._grids._Nx_left.value()
#                a2 = self._grids._Nx_left.value()*2
#                b = self._grids._Ny_left.value()
#                c = self._grids._Nz_left.value()
#                tem = 'potential_compass = "0 %d %d 0 %d 0 %d"\n'%(a1,a2,b,c)
#                tem += 'chargedensity_compass = "0 %d %d 0 %d 0 %d"\n'%(a1,a2,b,c)
#
#                with open(dir_name + '/input','w') as inc_file:
#                   inc_file.write('start_mode_NEGF="111"')
#                   inc_file.write(_mystate['input_grids_left3_lines'])
#                   inc_file.write(_mystate['input_mdscf_lines_NEGF'])
#                   inc_file.write(common_lines)
#                   inc_file.write(position)
#                   inc_file.write(tem)
#                   inc_file.write(_mystate['default_input_forON'])
#
#            except:
#               print "failed to save3a"
#            try:
#
#                conf1 = electrods[2]
#                conf2 = electrods[2]
#                conf3 = electrods[2]
#                position = position_orbital_3part(conf1,conf2,conf3,
#                                   num_orbital_dict, orbital_radius_dict)
#
#                dir_name = '3lead_lead2'
#                if not os.path.exists(dir_name): os.mkdir(dir_name)
#                a1 = self._grids._Nx_right.value()
#                a2 = self._grids._Nx_right.value()*2
#                b = self._grids._Ny_right.value()
#                c = self._grids._Nz_right.value()
#                tem = 'potential_compass = "0 %d %d 0 %d 0 %d"\n'%(a1,a2,b,c)
#                tem += 'chargedensity_compass = "0 %d %d 0 %d 0 %d"\n'%(a1,a2,b,c)
#
#                with open(dir_name + '/input','w') as inc_file:
#                   inc_file.write('start_mode_NEGF="111"')
#                   inc_file.write(_mystate['input_grids_right3_lines'])
#                   inc_file.write(_mystate['input_mdscf_lines_NEGF'])
#                   inc_file.write(common_lines)
#                   inc_file.write(position)
#                   inc_file.write(tem)
#                   inc_file.write(_mystate['default_input_forON'])
#
#            except:
#               print "failed to save3b"
#            try:
#                conf1 = electrods[1]
#                conf2 = electrods[0]
#                conf3 = electrods[2]
#                position = position_orbital_negf(conf1,conf2,conf3,
#                                   num_orbital_dict, orbital_radius_dict)
#
#                dir_name = 'bias_0.0'
#                if not os.path.exists(dir_name): os.mkdir(dir_name)
#                a1 = self._grids._Nx_left.value()
#                a2 = self._grids._Nx_negf-self._grids._Nx_right.value()
#                b = self._grids._Ny_right.value()
#                c = self._grids._Nz_right.value()
#                tem = 'potential_compass = "1 %d %d 0 %d 0 %d"\n'%(a1,a2,b,c)
#                tem += 'chargedensity_compass = "1 %d %d 0 %d 0 %d"\n'%(a1,a2,b,c)
#
#                with open(dir_name + '/input','w') as inc_file:
#                   inc_file.write('start_mode_NEGF="112"')
#                   inc_file.write(_mystate['input_grids_negf_lines'])
#                   inc_file.write(_mystate['input_mdscf_lines_NEGF'])
#                   inc_file.write(common_lines)
#                   inc_file.write(position)
#
#                   inc_file.write(tem)
#
#                   inc_file.write(_mystate['default_input_forON'])
#            except:
#               print "failed to save3c"
##            try:
##                        
##                write_out_LCR(self._io, self._grids, center_ctrl, configuration, num_orbital_dict)
##            except:
##               print "failed to save3d"
#            try:
#                write_out_others(configuration, self._rmg_para, self._grids, self._mdscf, num_orbital_dict)
#            except:
#               print "failed to save3e"
            try:
                write_out_jobfiles(configuration, self._setup, self._grids)
            except:
               print "failed to save job"
                    
        def state(self):
            """
               @return The state of all widgets in a state dictionary.
            """
            state = {}
            for widget in self._widgets:
#                print 's'
                state.update(widget.state())
#            print 'r'
            return state

        def title(self):
            """
               @return The title of the plugin.
            """
            return 'QW-RMG NCSU'
                
    except:
        showError()



    def selectdir(self):
        directory = QtGui.QFileDialog.getExistingDirectory(self)
        self.savedir.setText(directory)

        
def main():
    
    app = QtGui.QApplication(sys.argv)
    ex = RMG_GUI()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()

