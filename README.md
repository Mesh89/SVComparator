# SVComparator

Compare structural variations in a called VCF to a benchmark VCF, determining sensitivity and precision.

## Compiling the code

```
git clone https://github.com/Mesh89/SVComparator
cd SVComparator
./build_htslib.sh
cmake -DCMAKE_BUILD_TYPE=Release . && make
```

## Usage

Deletions are compared using compare-del, while insertions are compared using compare-ins.

```
compare-del $BENCH_VCF $CALLED_VCF $TRF_BED $REF $MAX_DIST [--report]
compare-ins $BENCH_VCF $CALLED_VCF $TRF_BED $REF $MAX_DIST [--ignore-seq --report]
```

Where 
- BENCH_VCF and CALLED_VCF are the benchmark and the called VCF files, respectively.
- TRF_BED are the TRF annotations for the used reference, in BED format
- REF is used reference in FASTA format
- MAX_DIST is the maximum distance between the SVs for them to be considered a match
- --report print precision and recall
- --ignore-seq does not compare the inserted sequences of two insertions


