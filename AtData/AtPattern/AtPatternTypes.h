#ifndef ATPATTERNTYPES_H
#define ATPATTERNTYPES_H

namespace AtPatterns {

/**
 * All implemented patterns. Can be created with the static factory AtPattern::CreatePattern(Type type)
 * @ingroup AtPattern
 */
enum class PatternType { kLine, kCircle2D }; //< Supported patterns

} // namespace AtPatterns
#endif //#ifndef ATPATTERNTYPES_H
