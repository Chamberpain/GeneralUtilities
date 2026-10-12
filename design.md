# GeneralUtilities design

## CEC support additions — 2026-10-12

Approved task scope: reusable path policy, native ROMS archives, and native rendering. Source code remains in this repository; CEC contains regional configuration and orchestration.

Data/Filepath/managed.py lazily constructs the existing handlers, permits only their out/tmp destinations, rejects data/store writes and escaping links, and fails when a configured root is unavailable. It never selects a fallback location.

Data/Download/native_roms_archive.py accepts a provider, geographic bounds, approved output and scratch directories, and bounded retry settings. WCOFS routing reuses the existing downloader; UCSC provides a collection. Selection retains nearest native C-grid windows, all sigma/interface levels, raw types, packed values, fill values and interpretation dependencies. Actual u/v time dimensions determine scheduling. Temporary computations stay in supplied scratch; output-side pending copies are atomically published. Immutable request identity, file checksums, advisory locks and manifest validation preserve existing valid outputs and support interrupted-write recovery. No regridding or time substitution occurs.

Data/Mapped/Utilities.py creates geographic axes and forwards native coordinate/value arrays directly to rendering methods. Project plotting implementations supply handler-derived output/cache directories.

Validation: 32 native archive tests, mapped-array preservation tests, the CEC integration suite and small real-source probes. Deployment was approved by Paul on 2026-10-12.

Final validation 2026-10-12: 476 tests and 10 subtests passed in combined CEC/GeneralUtilities integration. Real native probes passed; full exact-native mirror staged in approved scratch and removed after validation. Both missing WCOFS days are complete and bulk model saves continue. No staging, commit, merge or push.

Deployment scope: native ROMS writer, managed path policy, reusable mapping helpers, their tests and project documentation. The pre-existing depth_utilities.py edit and local AGENTS.md are excluded. User approved commit, merge and push; clean-tree verification precedes deployment.

Clean staged-tree verification: 60 checks passed with GeneralUtilities imported from an isolated ComputePathHandler snapshot and the pre-existing depth change excluded. Report: /Volumes/Compute/CEC/reports/out/deployment_generalutilities_tests_20261012.xml.
