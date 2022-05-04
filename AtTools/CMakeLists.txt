Set(LIBRARY_NAME AtTools)
Set(LINKDEF ${LIBRARY_NAME}LinkDef.h)

set(SRCS
  AtParsers.cxx
  AtELossManager.cxx

  AtCSVReader.cxx
  AtEulerTransformation.cxx
  AtSpaceChargeModel.cxx
  AtLineChargeModel.cxx
  AtVirtualTerminal.cxx
  AtKinematics.cxx

  AtFormat.cxx
  AtRandomSample.cxx
  )

Set(DEPENDENCIES
  ROOT::XMLParser
  ROOT::Core
  
  FairRoot::FairTools
  ATTPCROOT::AtData
)

generate_target_and_root_library(${LIBRARY_NAME}
  LINKDEF ${LINKDEF}
  SRCS ${SRCS}
  DEPS_PUBLIC ${DEPENDENCIES}
  )

