#ifndef ATMODELFACTORY_H
#define ATMODELFACTORY_H

#include "TString.h"
class AtTrackModel;

enum class AtModelType { kLINE };

class AtModelFactory {
public:
   static std::unique_ptr<AtTrackModel> CreateModel(AtModelType model);

   AtModelFactory() = delete;
   ~AtModelFactory() = delete;
};

#endif // #ifndef ATMODELFACTORY_H
