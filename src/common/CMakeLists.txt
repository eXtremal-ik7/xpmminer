include_directories(
  ${GMP_INCLUDE_DIRECTORY}
  ${OPENSSL_INCLUDE_DIR}
)

add_library(common STATIC
  CSieveOfEratosthenesL1Ext.cpp
  getblocktemplate.cpp
  primecoin.cpp
  system.cpp
)
