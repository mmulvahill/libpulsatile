# CI/CD Build Optimization Plan

## Goal
Reduce CI build time from **~18 minutes to 6-8 minutes** (60-65% reduction)

## Current Bottlenecks
1. **No caching** - R packages rebuild from source every time (5 min)
2. **RInside installation** - Always fails but wastes 2-4 min trying
3. **C++ standalone build** - Depends on RInside, usually skipped (wasted effort)
4. **System dependencies** - Reinstalled from scratch (2-3 min)
5. **apt-get update** - Unnecessary 1 min delay

## Implementation Strategy (4 Phases)

### Phase 1: Remove Waste (2-4 min savings)
**Effort:** 30 min | **Risk:** Very Low

**Changes:**
1. Delete RInside installation (ci.yml lines 53-71, coverage.yml lines 37-48)
   - R package uses `-DNORINSIDE` flag, doesn't need it
   - Installation fails 80% of time anyway

2. Delete C++ standalone build (ci.yml lines 73-95, coverage.yml lines 50-59)
   - Requires RInside (being removed)
   - R package compiles same C++ code via Rcpp

3. Remove `apt-get update` (ci.yml line 30, coverage.yml line 23)
   - GitHub Actions runners have recent packages
   - Saves ~1 min

**Expected Result:** 15-16 min builds after Phase 1

---

### Phase 2: R Package Caching (3-5 min savings)
**Effort:** 20 min | **Risk:** Low

**Add after "Setup R" step in both workflows:**

```yaml
- name: Cache R packages
  uses: actions/cache@v4
  with:
    path: ${{ env.R_LIBS_USER }}
    key: ${{ runner.os }}-r-${{ matrix.r-version }}-pkgs-v1-${{ hashFiles('R-package/DESCRIPTION') }}
    restore-keys: |
      ${{ runner.os }}-r-${{ matrix.r-version }}-pkgs-v1-
```

**How it works:**
- First run: Builds packages normally (~5 min), saves to cache
- Subsequent runs: Restores from cache (~20 sec)
- Invalidates when DESCRIPTION changes (new dependencies)
- Cache size: ~370 MB (well under 10GB limit)

**Expected Result:** 10-12 min builds on cache hit (70-80% of builds)

---

### Phase 3: System Dependency Caching (1-2 min savings)
**Effort:** 15 min | **Risk:** Low

**Add before system deps install (Ubuntu only):**

```yaml
- name: Cache system dependencies (Ubuntu)
  if: matrix.os == 'ubuntu-latest'
  uses: actions/cache@v4
  with:
    path: /var/cache/apt/archives
    key: ${{ runner.os }}-apt-v1-${{ hashFiles('.github/workflows/ci.yml') }}
    restore-keys: |
      ${{ runner.os }}-apt-v1-
```

**Expected Result:** 8-10 min builds

---

### Phase 4: Final Optimizations (0-1 min savings)
**Effort:** 15 min | **Risk:** Very Low

1. **Add parallel compilation:** Add `MAKEFLAGS: "-j2"` to env
2. **Streamline R install:** Single Rscript call with `upgrade='never'`

**Expected Result:** 7-9 min builds

---

## Implementation Order

1. **Phase 1** → Branch: `ci/remove-unused-steps`
   - Delete RInside, C++ build, apt-get update
   - Test: Verify all tests pass, ~15 min builds

2. **Phase 2** → Branch: `ci/add-r-package-caching`
   - Add R package cache
   - Test: First run ~15 min (cache build), second run ~10 min (cache hit)

3. **Phase 3** → Branch: `ci/add-apt-cache`
   - Add apt cache for Ubuntu
   - Test: Verify Ubuntu speedup, macOS unaffected

4. **Phase 4** → Branch: `ci/final-optimizations`
   - Add MAKEFLAGS and streamline install
   - Test: Final time ~7-9 min

---

## Critical Files

| File | Changes | Lines |
|------|---------|-------|
| `.github/workflows/ci.yml` | Remove RInside (53-71), C++ build (73-95), apt update (30)<br>Add R cache (after 25), apt cache (after 25)<br>Add MAKEFLAGS env, optimize install (44-51) | Major |
| `.github/workflows/coverage.yml` | Remove RInside (37-48), C++ build (50-59), apt update (23)<br>Add R cache (after 19), apt cache (after 19)<br>Optimize install (29-36) | Major |

---

## Rollback Strategy

**If caching breaks:**
1. Bump cache version: `pkgs-v1` → `pkgs-v2` (forces rebuild)
2. Remove cache step entirely if needed
3. Keep `ci.yml.backup` for emergency revert

**Detection:**
- "Package not found" errors → cache corruption
- Missing libraries → apt cache issue
- Both fixed by version bump or cache removal

---

## Success Metrics

**Primary:**
- Average build time: 7-9 min (target)
- Cache hit rate: >70%
- Test pass rate: 100% (unchanged)

**Secondary:**
- Time saved per week: ~50 hours across all builds
- Developer productivity: Faster PR feedback
- CI cost reduction: ~50% fewer compute minutes

---

## Testing Checklist

### Phase 1 (Remove Waste)
- [ ] Ubuntu build passes all tests
- [ ] macOS build passes all tests
- [ ] Build time <16 min
- [ ] No new warnings/errors

### Phase 2 (R Cache)
- [ ] First build creates cache (~15 min)
- [ ] Second build uses cache (~10 min)
- [ ] DESCRIPTION change rebuilds cache
- [ ] Cache restoration works correctly

### Phase 3 (Apt Cache)
- [ ] Ubuntu install time <1 min (vs 2-3 min)
- [ ] macOS unaffected
- [ ] All system libraries available

### Phase 4 (Final)
- [ ] Parallel compilation works (no race conditions)
- [ ] Build time 7-9 min average
- [ ] All optimizations stable

---

## Expected Timeline

- **Phase 1:** 30 min implementation + 20 min testing = 50 min
- **Phase 2:** 20 min implementation + 30 min testing = 50 min
- **Phase 3:** 15 min implementation + 20 min testing = 35 min
- **Phase 4:** 15 min implementation + 15 min testing = 30 min

**Total:** ~3 hours across 4 PRs

---

## Final Result

**Before:** 18 min builds, no caching, wasted effort on RInside
**After:** 7-9 min builds (cache hit), 15 min (cache miss), streamlined workflow

**Savings:** 50% average reduction, 60-65% on cache hit
