# QRSIDE 1.0.1

- Include the repository's deconvolution core in QRSIDE and implement spatial-radius neighbors internally; eliminate the conflicting SpatialDecon and SC.MEB installation requirements.
- Remove debug submatrix accesses that failed with embedding dimensions below five.
- Preserve matrix dimensions, cell-type names, spot alignment, and marker genes throughout preprocessing and model fitting.
- Validate raw counts, coordinates, marker formats, and model dimensions before expensive computation.
- Bound marker-stage iterations correctly, support a single spatial topic, and make initialization reproducible per call.
- Distribute a small source archive with complete synthetic and real MOB examples, and check the archive against the maintained source.
- Exercise installation, full model execution, and R package checks on Linux, macOS, and Windows.
