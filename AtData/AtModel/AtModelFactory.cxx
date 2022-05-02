#include "AtModelFactory.h"

#include "AtModelLine.h"

/**
 * @brief Create model.
 * Create and return a unique_ptr to a new model of the corresponding type
 * @param type Type of model to create
 * @return Pointer of new model
 */
std::unique_ptr<AtTrackModel> AtModelFactory::CreateModel(AtModelType type)
{
   switch (type) {
   case (AtModelType::kLINE): return std::make_unique<AtModelLine>();
   default: return nullptr;
   }
}
