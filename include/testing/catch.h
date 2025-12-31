//
// catch.h - Conditional Catch Framework Include
//
// When building standalone C++ tests (bin/tests), uses our bundled Catch v2.13.10
// When building R package tests (via testthat), uses testthat's vendor/catch.h
//
// This ensures test files work in both contexts with same syntax.
//

#ifndef GUARD_testing_catch_h
#define GUARD_testing_catch_h

#ifdef NORINSIDE
// R package build: use testthat's Catch framework
#include <testthat/vendor/catch.h>
#else
// Standalone C++ build: use our bundled Catch v2.13.10
#include <testing/catch_standalone.h>
#endif

#endif // GUARD_testing_catch_h
