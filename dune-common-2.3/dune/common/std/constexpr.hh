#ifndef DUNE_COMMON_STD_CONSTEXPR_HH
#define DUNE_COMMON_STD_CONSTEXPR_HH

#if HAVE_CONSTEXPR
#define DUNE_CONSTEXPR constexpr
#else // #if HAVE_CONSTEXPR
#define DUNE_CONSTEXPR
#endif // #else // #if HAVE_CONSTEXPR

#endif // #ifndef DUNE_COMMON_STD_CONSTEXPR_HH