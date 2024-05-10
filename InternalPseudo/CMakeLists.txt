add_library (InternalPseudo
InternalPseudo.cpp
InternalPseudo_uspp.cpp
InternalPseudo_sg15.cpp
InternalPseudo_ncpp_stringent.cpp
InternalPseudo_ncpp_standard.cpp
)
include_directories("${RMG_SOURCE_DIR}/InternalPseudo/Data")

