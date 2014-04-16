include_directories(
  ${GMP_INCLUDE_DIRECTORY}
  ${OPENSSL_INCLUDE_DIRECTORY}
  ${JANSSON_INCLUDE_DIRECTORY}
  ${CURL_INCLUDE_DIRECTORY}
)

add_executable(xpmcpuminer
  xpmcpuminer.cpp
  rippedFromHp.cpp
) 

target_link_libraries(xpmcpuminer
  common
  blkmaker
  ${GMP_LIBRARY}
  ${GMPXX_LIBRARY}
  ${SSL_LIBRARY}
  ${CRYPTO_LIBRARY}
  ${CURL_LIBRARY}
  ${JANSSON_LIBRARY}
  ${NCURSES_LIBRARY}
)

if (OS_LINUX)
  target_link_libraries(xpmcpuminer pthread)
endif()

install(TARGETS xpmcpuminer DESTINATION bin)
