add_library(XClibs
    dft_setting_routines.F90
    dft_setting_params.F90
    qe_constants.F90
    qe_dft_list.F90
    qe_drivers_d_gga.F90
    qe_drivers_d_lda_lsda.F90
    qe_drivers_gga.F90
    qe_drivers_lda_lsda.F90
    qe_drivers_mgga.F90
    qe_funct_corr_gga.F90
    qe_funct_corr_lda_lsda.F90
    qe_funct_exch_gga.F90
    qe_funct_exch_lda_lsda.F90
    qe_funct_mgga.F90
    qe_kind.F90
    xclib_error.F90
    xclib_utils_and_para.F90
#    xc_beef_interface.F90
    xc_lib.F90
    xc_wrapper_d_gga.F90
    xc_wrapper_d_lda_lsda.F90
    xc_wrapper_gga.F90
    xc_wrapper_lda_lsda.F90
    xc_wrapper_mgga.F90)
    
include_directories("${RMG_BINARY_DIR}/XC/")
include_directories("${RMG_BINARY_DIR}/XClibs/")
