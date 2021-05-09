include_directories(
  ${GMP_INCLUDE_DIRECTORY}
  ${OPENSSL_INCLUDE_DIR}
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
  ${OPENSSL_CRYPTO_LIBRARY}
  ${OPENSSL_SSL_LIBRARY}
  ${CURL_LIBRARY}
  ${JANSSON_LIBRARY}
)

if (OS_LINUX)
  target_link_libraries(xpmcpuminer pthread dl)
elseif (OS_WINDOWS)
  set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++ -static")
  target_link_libraries(xpmcpuminer
    ${PTHREAD_LIBRARY}
    ${Z_LIBRARY}
    ws2_32
    wldap32
  )
endif()

install(TARGETS xpmcpuminer DESTINATION bin)
