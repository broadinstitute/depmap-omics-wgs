# Annotate mutations postprocessing

This module loads the DuckDB database generated by [vcf-to-duckdb](https://github.com/broadinstitute/vcf-to-duckdb) and exports a single Parquet file representing a MAF file. This file contains only the filtered structural variants and only the columns needed by the [postprocessing function](https://github.com/broadinstitute/depmap_omics/blob/master/depmapomics/mutations.py#L440) in depmap-omics.
