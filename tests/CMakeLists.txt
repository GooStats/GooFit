file(GLOB_RECURSE test_cases *.cc)
foreach(case_file ${test_cases})
get_filename_component( case_name ${case_file} NAME_WE )
set (case_name test_${case_name})
add_executable(${case_name} ${case_file})
  target_link_libraries(${case_name} goofit_lib
      gtest gtest_main)
  add_test(NAME ${case_name}
      COMMAND $<TARGET_FILE:${case_name}>
      WORKING_DIRECTORY
      ${PROJECT_BINARY_DIR})
endforeach()
