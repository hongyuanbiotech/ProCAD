# ProCAD

ProCAD is an R script to normalized coverage data from next-generation-sequencing data. 


## Usage

```bash
R --no-save zscore.R 'the folder contains raw coverage data'.
```

The coverage data should be pre-formatted as following
chrom start_position  raw_coverage_in_doubles
chrX    155200000       5.10025
.....

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
