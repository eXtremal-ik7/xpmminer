include_directories(${GMP_INCLUDE_DIRECTORY})

add_library(common STATIC
  CSieveOfEratosthenesL1Ext.cpp
  getblocktemplate.cpp
  primecoin.cpp
  system.cpp
)
