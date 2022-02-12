add_library(GooFitTestFit fit.cc)
target_link_libraries(GooFitTestFit PUBLIC goofit_lib)
add_executable(test_GooFit expFit.cc)
target_link_libraries(test_GooFit PRIVATE GooFitTestFit gtest_main)
target_include_directories(test_GooFit PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})


add_test(NAME test_GooFit
  COMMAND $<TARGET_FILE:test_GooFit>)

add_library(GooFitTestFit_NLL fit.cc)
target_link_libraries(GooFitTestFit_NLL PUBLIC goofit_lib_NLL)
add_executable(test_GooFitNLL NLLtest.cc)
target_link_libraries(test_GooFitNLL PRIVATE GooFitTestFit_NLL gtest_main ROOT::RIO)
target_include_directories(test_GooFitNLL PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
file(COPY NLL_CHECK_reference.root DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_GooFitNLL
      COMMAND $<TARGET_FILE:test_GooFitNLL>)