goofit_add_library(PDFs STATIC
  PdfBase.cu
  GooPdf.cu
  FitControl.cu
  SumPdf.cu
  ExpPdf.cu
)

target_link_libraries(PDFs PUBLIC goofit_base MyThrust)
