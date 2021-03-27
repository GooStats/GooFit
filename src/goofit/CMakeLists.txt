add_library(goofit_base
	BinnedDataSet.cc
	DataSet.cc
	FitControl.cc
	FitManager.cc
	PdfBase.cc
	UnbinnedDataSet.cc
	Variable.cc
	)

target_link_libraries(goofit_base PUBLIC ROOT::Core ROOT::Minuit)
target_compile_features(goofit_base PUBLIC cxx_std_14)
