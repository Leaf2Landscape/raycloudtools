set(SOURCES
  rayvoxel.cpp
)

ras_add_executable(rayvoxel
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
