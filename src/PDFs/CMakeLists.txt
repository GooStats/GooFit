goofit_add_library(OtherPDFs GooPdf.cu PdfBase.cu FitControl.cu ExpPdf.cu)
target_link_libraries(OtherPDFs PUBLIC goofit_base MyThrust)

goofit_add_library(SumPdf SumPdf.cu)
target_link_libraries(SumPdf PUBLIC goofit_base MyThrust)

add_library(PDFs INTERFACE)
target_link_libraries(PDFs INTERFACE SumPdf OtherPDFs)

goofit_add_library(SumPdf_NLL SumPdf.cu)
target_compile_options(SumPdf_NLL PUBLIC "-DNLL_CHECK=ON")
target_link_libraries(SumPdf_NLL PUBLIC goofit_base MyThrust GooStatsNLLCheck)

add_library(PDFs_NLL INTERFACE)
target_link_libraries(PDFs_NLL INTERFACE SumPdf_NLL OtherPDFs)
