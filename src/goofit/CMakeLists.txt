add_library(goofit_base
#	BinnedDataSet.cc
#	DataSet.cc
	FitControl.cc
	FitManager.cc
	PdfBase.cc
#	UnbinnedDataSet.cc
#	Variable.cc
	)

find_package(ROOT CONFIG REQUIRED COMPONENTS Minuit TreePlayer)
target_link_libraries(goofit_base PUBLIC ROOT::Core)
target_compile_features(goofit_base PUBLIC cxx_std_14)
message(STATUS ${PROJECT_SOURCE_DIR}/include)
