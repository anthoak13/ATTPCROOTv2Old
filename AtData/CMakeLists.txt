# Create a library called "libAtReco" which includes the source files given in the array .

Set(LIBRARY_NAME AtData)
Set(LINKDEF AtDataLinkDef.h)


set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
)

set(LINK_DIRECTORIES
${FAIRROOT_LIBRARY_DIR}
)

Set(DEPENDENCIES
  Base
  FairTools

  #AtReconstruction

  ROOT::Spectrum
  ROOT::Core
  ROOT::GenVector
  )

set(LINK_DIRECTORIES
  ${FAIRROOT_LIBRARY_DIR}
  )

set(SRCS
  # Add all the source files below this line. Those must have cxx for their extension.
  AtHit.cxx
  AtPad.cxx
  AtRawEvent.cxx
  AtEvent.cxx
  AtHitCluster.cxx
  AtTrack.cxx
  AtTrackingEventAna.cxx
  AtProtoQuadrant.cxx
  AtProtoEvent.cxx
  AtProtoEventAna.cxx
  AtPatternEvent.cxx
  )

generate_target_root_library(${LIBRARY_NAME}
  LINKDEF ${LINKDEF}
  SRCS ${SRCS}
  INCLUDE_DIR ${INCLUDE_DIRECTORIES}
  LIBRARY_DIR ${LINK_DIRECTORIES}
  DEPS_PUBLIC ${DEPENDENCIES}
  )
