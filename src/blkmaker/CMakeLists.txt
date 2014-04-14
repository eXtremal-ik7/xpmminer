include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${JANSSON_INCLUDE_DIRECTORY}
)

set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} "-std=c99")

add_library(blkmaker STATIC
  base58.c
  blkmaker.c
  blkmaker_jansson.c
  blktemplate.c
  hex.c
)
