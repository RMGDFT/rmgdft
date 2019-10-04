
from PyQt5 import QtGui, QtCore

from distutils.sysconfig import get_python_lib

## presumably need to import pickle module twice to make it work properly:
import json
import codecs

class default_input_para():
    """
       Widget for default input para
    """


    def state(self):
        """
           @return A dictionary containing the widget state.
        """



        default_input_forON  = '\n#  **** default parameters for order-n input  ***** \n\n'
        default_input_forON += 'kpoints_per_processor="1"\n'
        default_input_forON += 'potential_grid_refinement="2"\n'
        default_input_forON += 'beta_grid_refinement="4"\n'
        default_input_forON += 'calculation_mode="Quench Electrons"\n'
        default_input_forON +='number_of_kpoints="1"\n'
        default_input_forON += 'kpoints= "0.0  0.0  0.0   1.0 "\n'
        default_input_forON += 'do_spin_polarized="false"\n'
        default_input_forON += 'energy_cutoff_parameter="1.75"\n'
        default_input_forON += 'kohn_sham_pre_smoothing="3"\n'
        default_input_forON += 'kohn_sham_post_smoothing="3"\n'
        default_input_forON += 'poisson_pre_smoothing="2"\n'
        default_input_forON += 'poisson_post_smoothing="2"\n'
        default_input_forON += 'do_write_waves_to_file="true"\n'
        default_input_forON += 'fine_grid_non_local_pp="4"\n'
        default_input_forON += 'mg_method="Pulay"\n'
        default_input_forON += 'mg_steps="2"\n'
        default_input_forON += 'do_movable_orbital_centers="false"\n'
        default_input_forON += 'movable_orbital_centers_steps="1"\n'

        default_input_for_oneatom  = '\n#  **** default parameters for one atom input  ***** \n\n'
        default_input_for_oneatom +='num_processor ="1"\n'
        default_input_for_oneatom +='processor_grid="1 1 1"\n'
        default_input_for_oneatom +='kpoints_per_processor="1"\n'
        default_input_for_oneatom +='Hamiltonia_processor_grid="1 1"\n'
        default_input_for_oneatom +='start_mode="Random Start"\n'
        default_input_for_oneatom +='number_of_atoms="1"\n'
        default_input_for_oneatom +='number_of_species="1"\n'

        state={ 'default_input_for_oneatom': default_input_for_oneatom,
                'default_input_forON': default_input_forON
                 }
        return state

    ######### end of state(self):
