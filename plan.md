# GeneralUtilities plan

## CEC native-data and file-layout support — 2026-10-12

Status: approved by Paul with the CEC tasks; build/test completed, deployment approved by Paul on 2026-10-12. Branch: cec/native-layout. Preserve the pre-existing depth_utilities.py edit and other unrelated work.

1. Add lazy output/scratch subclasses of the existing FilePathHandler and ComputePathHandler without implicit shared-data/store creation or alternate destinations.
2. Add a provider-based native ROMS archive writer under Data/Download. Preserve raw C-grid fields, source metadata and sigma layers; use nearest native index windows and exact source time selection.
3. Add small native-array rendering adapters under Data/Mapped for project plotting scripts.
4. Test stored-value/metadata preservation, geometry, lock/resume/conflict handling, cross-filesystem publication, path constraints, and unaltered rendering inputs. Tests live in unit_tests.
5. Confirm with Paul before deployment commit, merge, or push.

The regional acquisition schedules, real probes, output migration and operational status are documented in the CEC project plan/design and handler-root reports. No shared HYCOM download behavior is changed.

Final validation 2026-10-12: 476 tests and 10 subtests passed in combined CEC/GeneralUtilities integration. Real native probes passed; full exact-native mirror staged in approved scratch and removed after validation. Both missing WCOFS days are complete and bulk model saves continue. No staging, commit, merge or push.

Deployment scope: native ROMS writer, managed path policy, reusable mapping helpers, their tests and project documentation. The pre-existing depth_utilities.py edit and local AGENTS.md are excluded. User approved commit, merge and push; clean-tree verification precedes deployment.

Clean staged-tree verification: 60 checks passed with GeneralUtilities imported from an isolated ComputePathHandler snapshot and the pre-existing depth change excluded. Report: /Volumes/Compute/CEC/reports/out/deployment_generalutilities_tests_20261012.xml.
