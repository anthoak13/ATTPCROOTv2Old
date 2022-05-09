#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace AtTools;

#pragma link C++ class AtTools::AtELossManager + ;
#pragma link C++ class AtTools::AtParsers + ;
#pragma link C++ class AtEulerTransformation + ;

#pragma link C++ class AtSpaceChargeModel + ;
#pragma link C++ class AtLineChargeModel + ;

#pragma link C++ class AtTools::AtKinematics + ;
#pragma link C++ class AtTools::AtVirtualTerminal + ;

#pragma link C++ class AtTools::AtSample - !;
#pragma link C++ class AtTools::AtIndependentSample - !;
#pragma link C++ class AtTools::AtUniform - !;
#pragma link C++ class AtTools::AtChargeWeighted - !;
#pragma link C++ class AtTools::AtGaussian - !;
#pragma link C++ class AtTools::AtWeightedGaussian - !;

#endif
