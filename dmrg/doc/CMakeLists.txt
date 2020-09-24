include(UseLATEX)
add_latex_document(qcmaquis_manual.tex BIBFILES qcmaquis_manual.bib)
install(FILES ${CMAKE_BINARY_DIR}/doc/qcmaquis_manual.pdf COMPONENT doc DESTINATION ${DOC_INSTALL_PATH})
