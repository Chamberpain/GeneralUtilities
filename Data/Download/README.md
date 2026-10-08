# Model downloaders

Ocean model transport, subset validation, retry logic and archive publication are
maintained here. HyperNav, CEC and SeatrecPlots retain their regional configuration,
plotting and simulation adapters. Their original import paths and CLI scripts
continue to work through the shared implementations.

| Module | Purpose |
| --- | --- |
| `hycom_download` | `HYCOMDownloader`: ESPC source discovery, regional metadata, validated legacy pickle chunks and manifests; configure regional bounds, depth and a file handler in a subclass. |
| `hycom_forecast_download` | `return_tds_dataset` and `NCSS_dataset`: forecast U/V or T/S subsets, native depths and forecast-run precedence. |
| `hycom_october_download` | `archive_sources`, `hycom_archive_times` and `download_hycom_october`: historical source routing, adaptive requests and resumable October subsets. |
| `production_hycom_download` | Legacy production `Source`, `SourceNOAA`, `SourceFNMOCBest` and `SourceFNMOCRun` providers. |
| `copernicus_global_download` | `CopernicusDownloader`: lazy SDK subsets, canonical U/V/W archives and existing HyperNav pickle downloads. Regional classes supply bounds and location. |
| `copernicus_download` | `download_copernicus_data(start, end, bounds=(west, east, south, north), ...)`: validated regional SDK subset; an optional consumer validator runs before output publication. |
| `wcofs_download` | `WCOFSDownloader`: time-correct source URLs, physical-grid normalization and archive downloading. |
| `wcofs_roms` | Native ROMS destaggering, velocity rotation and depth/geographic interpolation. |
| `current_archive` | `download_record(fetch_chunk, ...)`: concurrent chunk scheduling, request identity, file locking, validation and atomic manifests. |
| `cec_hycom_download` | Existing archive U/V reader and missing-file repair, including template/grid validation. |
| `hycom_redownload`, `hycom_redownload_resilient` | Prepare/worker/merge repair CLIs; resilient requests add coordinated throttling and single-time fallback. |
| `legacy_model_download` | Shared legacy PyDAP/Motu pickle chunk loops. |

The forecast, historical and canonical archive interfaces deliberately retain
separate depth/time selection rules and file formats. Existing filenames,
manifest schemas, CLI defaults and explicit `--module-dir` overrides are
preserved. The unchanged default WCOFS ROMS processor retains its previous
provenance identity so prior archives can resume.

SDK credentials remain caller- or SDK-managed. Importing the shared engines does
not authenticate or download data. Service-specific dependencies are required
when their download methods are invoked; native NetCDF subset handling uses
`netCDF4`. Legacy cache adapters use `FilePathHandler` lazily and initialize only
`tmp` and `out`; use of `store` requires explicit authorization and a separately
configured handler. Explicit output directories are still supported.

Offline regressions are in the project-level `unit_tests/` folder. Existing HyperNav archive/normalization
and CEC repair regressions also exercise the compatibility entry points. Live
service requests are not part of these tests.
