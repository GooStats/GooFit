add_library(goofit_base
	BinnedDataSet.cc
	DataSet.cc
	FitControl.cc
	FitManager.cc
	PdfBase.cc
	UnbinnedDataSet.cc
	Variable.cc
	)

target_link_libraries(goofit_base PUBLIC MyThrust ROOT::Core ROOT::Minuit)
target_include_directories(goofit_base PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_compile_features(goofit_base PUBLIC cxx_std_14)
