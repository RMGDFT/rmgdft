class default_input_para:
    """
    Widget for default input para
    """

    def state(self):
        """
        @return A dictionary containing the widget state.
        """

        default_input_forON = (
            "\n#  **** default parameters for order-n input  ***** \n\n"
            'kpoints_per_processor="1"\n'
            'potential_grid_refinement="2"\n'
            'beta_grid_refinement="4"\n'
            'calculation_mode="Quench Electrons"\n'
            'number_of_kpoints="1"\n'
            'kpoints= "0.0  0.0  0.0   1.0 "\n'
            'do_spin_polarized="false"\n'
            'energy_cutoff_parameter="1.75"\n'
            'kohn_sham_pre_smoothing="3"\n'
            'kohn_sham_post_smoothing="3"\n'
            'poisson_pre_smoothing="2"\n'
            'poisson_post_smoothing="2"\n'
            'do_write_waves_to_file="true"\n'
            'fine_grid_non_local_pp="4"\n'
            'mg_method="Pulay"\n'
            'mg_steps="2"\n'
            'do_movable_orbital_centers="false"\n'
            'movable_orbital_centers_steps="1"\n'
        )

        default_input_for_oneatom = (
            "\n#  **** default parameters for one atom input  ***** \n\n"
            'num_processor ="1"\n'
            'processor_grid="1 1 1"\n'
            'kpoints_per_processor="1"\n'
            'Hamiltonia_processor_grid="1 1"\n'
            'start_mode="Random Start"\n'
            'number_of_atoms="1"\n'
            'number_of_species="1"\n'
        )
        state = {
            "default_input_for_oneatom": default_input_for_oneatom,
            "default_input_forON": default_input_forON,
        }
        return state
