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

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_INSTALL_LIBDIR .)
ROOT_GENERATE_DICTIONARY(G__GooStatsNLLCheck ${PROJECT_SOURCE_DIR}/include/GooStatsNLLCheck.h
  LINKDEF ${PROJECT_SOURCE_DIR}/include/GooStatsNLLCheckLinkDef.h)
add_library(GooStatsNLLCheck SHARED GooStatsNLLCheck.cc G__GooStatsNLLCheck)
target_include_directories(GooStatsNLLCheck PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_link_libraries(GooStatsNLLCheck PUBLIC ROOT::Core ROOT::RIO)
install(TARGETS GooStatsNLLCheck DESTINATION .)
