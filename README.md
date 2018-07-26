kseq++
======
kseq++ is a C++ re-implementation of [kseq](https://github.com/attractivechaos/klib/blob/master/kseq.h)
by [Heng Li](https://github.com/lh3). The goal for re-implementation of `kseq` is
providing better API and resource management while preserving the flexibility,
and performance. Like original kseq, this parser is based on generic stream
buffer and works with different file types. However, instead of using C macros,
it uses C++ templates. The RAII-class `KStream` is the main class which is
constructed by `make_kstream` function. It gets the file object (can be of any
type) and its corresponding read function. In contrast with kseq, there is no
need to specify the types, since they are inferred by compiler. Each record will
be stored in a `KSeq` object.

It inherits all features from kseq (quoting from kseq homepage):
> - Parse both FASTA and FASTQ format, and even a mixture of FASTA and FASTQ records in one file.
> - Seamlessly adapt to gzipped compressed file when used with zlib.
> - Support multi-line FASTQ.
> - Work on a stream with an internal stream buffer. 

Reading a sequence file
-----------------------

```c++
#include <iostream>
#include <zlib>
#include "kseq++.h"

int main(int argc, char* argv[])
{
  KSeq record;
  gzFile fp = gzopen(filename, "r");
  auto ks = make_kstream(fp, gzread);
  while (ks >> record) {
    std::cout << record.name << std::endl;
    if (!record.comment.empty()) std::cout << record.comment << std::endl;
    std::cout << record.seq << std::endl;
    if (!record.qual.empty()) std::cout << record.qual << std::endl;
  }
}
```

Benchmark
---------
### Datasets
For this benchmark, I re-used sequence files from SeqKit benchmark:
[seqkit-benchmark-data.tar.gz](http://app.shenwei.me/data/seqkit/seqkit-benchmark-data.tar.gz)

| file         | format | type |  num_seqs |       sum_len | min_len |      avg_len |     max_len |
| :----------- | :----- | :--- | --------: | ------------: | ------: | -----------: | ----------: |
| dataset_A.fa | FASTA  | DNA  |    67,748 | 2,807,643,808 |      56 |     41,442.5 |   5,976,145 |
| dataset_B.fa | FASTA  | DNA  |       194 | 3,099,750,718 |     970 | 15,978,096.5 | 248,956,422 |
| dataset_C.fq | FASTQ  | DNA  | 9,186,045 |   918,604,500 |     100 |          100 |         100 |

### Platform

- CPU: Intel&reg; Core&reg; i5-4200U CPU @ 1.60GHz, two cores, 4 threads
- RAM: DDR3 1600 MHz, 8192 MB
- SSD: Samsung SSD 850 EVO 500GB, SATA-3
- OS: Debian GNU/Linux 9.4 (stretch), Linux 4.9.0-6-amd64
- Compiler: GCC 6.3.0, compiled with optimisation level 3 (`-O3`)

### Result

| file         | kseq++ |   kseq |  SeqAn |
| :----------- | -----: | -----: | -----: |
| dataset_A.fa | 7.45 s | 6.22 s | 9.96 s |
| dataset_B.fa | 6.78 s | 5.71 s | 19.9 s |
| dataset_C.fq | 5.73 s | 4.21 s | 8.66 s |
