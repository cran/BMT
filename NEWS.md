# BMT 0.1.3

## Minor changes

- Corrected spelling errors in documentation files.
- Corrected the `Version` field in the DESCRIPTION file to follow the required `X.Y.Z` format.
- Added the `Authors@R` field to DESCRIPTION and removed deprecated `Author` and `Maintainer` fields.
- Set `Encoding: UTF-8` in DESCRIPTION as required by `roxygen2`.
- Added missing `@export` tags to functions intended to be public.
- Suppressed export of internal S3 method `startargdefault` by omitting the `@export` tag.
- Reorganized roxygen blocks to avoid multiple `@rdname` tags in the same block.
- Fixed invalid URLs and replaced direct `arXiv:` references with proper DOIs.
