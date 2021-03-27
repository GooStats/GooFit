goofit_add_library(PDFs STATIC
  PdfBase.cu
)

target_link_libraries(PDFs goofit_base MyThrust)
