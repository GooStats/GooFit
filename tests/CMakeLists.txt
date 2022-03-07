include(GoogleTest)

add_library(GooFitTestFit fit.cc)
add_sanitizers(GooFitTestFit)
target_link_libraries(GooFitTestFit PUBLIC goofit_lib)

add_executable(test_GooFit expFit.cc)
add_sanitizers(test_GooFit)
target_link_libraries(GooFitTestFit PUBLIC goofit_lib)
target_link_libraries(test_GooFit PRIVATE GooFitTestFit gtest_main)
target_include_directories(test_GooFit PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
gtest_discover_tests(test_GooFit)

add_library(GooFitTestFit_NLL fit.cc)
add_sanitizers(GooFitTestFit_NLL)
target_link_libraries(GooFitTestFit_NLL PUBLIC goofit_lib_NLL)
add_executable(test_GooFitNLL NLLtest.cc)
add_sanitizers(test_GooFitNLL)
target_link_libraries(test_GooFitNLL PRIVATE GooFitTestFit_NLL gtest_main
        ROOT::RIO)
target_include_directories(test_GooFitNLL PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
file(COPY NLL_CHECK_reference.root DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
gtest_discover_tests(test_GooFitNLL)
