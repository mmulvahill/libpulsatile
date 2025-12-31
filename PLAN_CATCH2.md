# Plan: Evaluate and Merge Changes from develop Branch

## Analysis Summary

### Current State
- **Current branch**: `fix/deprecated-tidyverse-functions` (6 commits ahead of master)
- **develop branch**: 6 commits ahead of master
- **Uncommitted changes**: Modified CI/CD workflows with `Sys.getenv()` approach

### Key Finding: Current Branch is MORE Modern Than develop!

The R code in `fix/deprecated-tidyverse-functions` is actually **more up-to-date** than develop:

| Feature | develop (OLD) | Current Branch (NEW) |
|---------|--------------|---------------------|
| **tidyr** | `gather_()`, `spread()` | `pivot_longer()`, `pivot_wider()` |
| **ggplot2** | `aes_string()`, `size`, `..density..` | `aes()` with `.data$`, `linewidth`, `after_stat()` |
| **tibble** | `as_data_frame()` | `as_tibble()` |
| **Class checking** | `class(x) == "type"` | `inherits(x, "type")` |
| **R CMD check** | Missing globalVariables | Added `utils::globalVariables("density")` |

### What to Keep from develop Branch

#### 1. ✅ Catch2 Upgrade (CRITICAL - Must Keep)
**File**: `include/testing/catch.h`
- Upgrades from v2.0.1 (2017) → v2.13.10 (2022)
- **Critical fix**: SIGSTKSZ compatibility with modern glibc
- Prevents build failures on modern systems

#### 2. ✅ Test Improvements (Keep)
**File**: `tests/mh_tests.cpp`
- Adds `TestConfig` namespace with calibrated constants
- Improves test determinism and reduces flakiness
- Better floating-point comparison with `Approx()`

#### 3. ❌ Tidyverse Changes (DO NOT KEEP)
**Files**: `R-package/R/*.R`
- develop has OLDER, deprecated functions
- Current branch is already modernized

#### 4. ⚠️ CI/CD Changes (Needs Decision)
**Files**: `.github/workflows/ci.yml`, `.github/workflows/coverage.yml`

**develop approach**:
```bash
Rscript -e ".libPaths('$R_LIBS_USER'); install.packages(...)"
```
Uses shell variable expansion

**Current branch uncommitted approach**:
```bash
Rscript -e ".libPaths(Sys.getenv('R_LIBS_USER')); install.packages(...)"
```
Uses R's environment variable access

Both fix the same issue but with different methods.

## Questions for User

Need clarification on the following before proceeding:

1. **CI/CD Approach**: Which approach do you prefer?
   - A) Keep develop's shell variable expansion approach
   - B) Keep your uncommitted `Sys.getenv()` approach
   - C) Test both and see which works better

2. **Branch Strategy**: How should we merge?
   - A) Cherry-pick only Catch2 + test improvements from develop → current branch
   - B) Merge develop into current branch, then revert R code to current versions
   - C) Start fresh: merge master → current branch, then apply both sets of changes

3. **Other changes from develop**: Do you want to keep?
   - CLAUDE.md updates
   - README.md formatting cleanup
   - RoxygenNote version bump (6.1.1 → 7.3.1)
   - NAMESPACE changes (import updates)

## User Decisions

✅ **CI/CD**: Make NO changes - another Claude instance is managing this
✅ **Merge strategy**: Cherry-pick specific commits from develop
✅ **Metadata files**: Keep RoxygenNote bump, NAMESPACE changes, README cleanup

## Final Analysis: Current Branch is MORE Modern!

After detailed comparison, **the current branch is superior to develop** in almost every way:

| Component | develop | Current Branch | Winner |
|-----------|---------|---------------|--------|
| **Catch2** | v2.13.10 (2022) ✅ | v2.0.1 (2017) | develop |
| **Tests** | Disabled/reverted | Improved + deterministic ✅ | current |
| **R tidyverse** | OLD deprecated functions | NEW modern functions ✅ | current |
| **RoxygenNote** | 6.1.1 (OLD) | 7.3.1 (NEW) ✅ | current |
| **NAMESPACE** | OLD imports | NEW imports ✅ | current |
| **CLAUDE.md** | Better documentation ✅ | Basic | develop |
| **README.md** | Line wrapping cleanup ✅ | 78 char lines | develop |

## Implementation Plan

### Step 1: Cherry-pick Catch2 Upgrade from develop

**Commit**: `95f87fd` - "updated to newest version of catch2 (2.13.10)"

**Files modified by this commit**:
- `include/testing/catch.h` - Catch v2.0.1 → v2.13.10
- `CLAUDE.md` - Initial documentation about Catch2

**Action**:
- Cherry-pick this commit
- Resolve any conflicts in `CLAUDE.md` (prefer keeping current content + adding Catch2 details)

### Step 2: Merge CLAUDE.md Improvements from develop

Since `CLAUDE.md` is touched by multiple commits, manually merge the improvements:

**Changes to add from develop**:
```markdown
## Dependencies
- Enhanced R package list with categorization (Required/Development/Optional)
- Updated install command with full package list
- Explicit mention of RInside as optional

## Code Style Guidelines
- Update: "Tests: Catch2 v2.13.10 framework for C++"

## Common Issues
- Add: "ARM64 Mac compatibility: Upgraded to Catch2 v2.13.10 with native Apple Silicon support"
- Add: "Ubuntu/glibc compatibility: Modern Catch2 resolves SIGSTKSZ compilation errors"
- Add: "RcppArmadillo warnings: Added `-DARMA_USE_CURRENT` flag to suppress C++11 deprecation warnings"
- Add: "Testing framework: Upgraded from Catch v2.0.1 (2017) to Catch2 v2.13.10 (2022)"
```

### Step 3: Apply README.md Line Wrapping

**Change**: Line 16 wraps at 80 characters:
```markdown
# Before (current):
The package is currently in development, with only the single subject model currently functioning. Feel free to take a look around and test it out.

# After (develop):
The package is currently in development, with only the single subject model currently functioning. Feel free to
take a look around and test it out.
```

### Step 4: Keep Everything Else from Current Branch

**DO NOT take from develop**:
- ❌ `tests/mh_tests.cpp` - develop has disabled tests; current has improved tests
- ❌ `R-package/R/*.R` files - develop has OLD tidyverse; current has NEW
- ❌ `R-package/DESCRIPTION` - develop has OLD RoxygenNote 6.1.1; current has NEW 7.3.1
- ❌ `R-package/NAMESPACE` - develop has OLD imports; current has NEW
- ❌ `.github/workflows/*` - per user request, don't touch CI/CD
- ❌ `Makefile` - current branch is fine

## Execution Steps

1. **Cherry-pick Catch2 commit**:
   ```bash
   git cherry-pick 95f87fd
   ```

2. **Resolve CLAUDE.md conflict** (if any):
   - Keep current branch's CI/CD documentation
   - Add develop's enhanced dependency documentation
   - Add develop's Catch2 and platform compatibility notes

3. **Manually update CLAUDE.md** with improvements from develop:
   - Enhanced Dependencies section
   - Catch2 version specification
   - Platform compatibility notes

4. **Apply README.md line wrapping**:
   - Wrap line 16 at ~80 characters

5. **Verify no CI/CD files changed**:
   ```bash
   git status
   # Should NOT show .github/workflows/ changes
   ```

6. **Test the changes**:
   ```bash
   make clean && make
   bin/tests
   cd R-package && Rscript -e 'devtools::check()'
   ```

## Files to Modify

### From develop (via cherry-pick):
- `include/testing/catch.h` (automatic via cherry-pick)

### Manual edits:
- `CLAUDE.md` (merge improvements from develop)
- `README.md` (line wrapping)

### Keep unchanged from current branch:
- `tests/mh_tests.cpp` ✅
- `R-package/R/*.R` ✅
- `R-package/DESCRIPTION` ✅
- `R-package/NAMESPACE` ✅
- `.github/workflows/*` ✅ (per user request)

## Summary

This plan **selectively takes only the valuable parts** from develop (Catch2 upgrade + documentation) while **preserving all the modern improvements** in the current branch (modern tidyverse, better tests, newer roxygen).
