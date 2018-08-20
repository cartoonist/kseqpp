kseq++
======
kseq++ is a C++11 re-implementation of [kseq](https://github.com/attractivechaos/klib/blob/master/kseq.h)
by [Heng Li](https://github.com/lh3). The goal for re-implementation of `kseq` is
providing better API and resource management while preserving its flexibility
and performance. Like original kseq, this parser is based on generic stream
buffer and works with different file types. However, instead of using C macros,
it uses C++ templates. The RAII-style class `KStream` is the main class which
can be constructed by `make_kstream` function series or by calling its
constructor directly (C++17). It gets the file object/pointer (can be of any
type), its corresponding read/write function, and opening mode (`mode::in` or
`mode::out`).  In contrast with kseq, there is no need to specify the types,
since they are inferred by compiler. Each record will be stored in a `KSeq`
object.

It inherits all features from kseq (quoting from kseq homepage):
> - Parse both FASTA and FASTQ format, and even a mixture of FASTA and FASTQ records in one file.
> - Seamlessly adapt to gzipped compressed file when used with zlib.
> - Support multi-line FASTQ.
> - Work on a stream with an internal stream buffer.

while additionally provides:
- simpler and more readable API
- RAII-style memory management

The library also comes with FASTA/Q writer. Like reading, it can write mixed
multi-line FASTA and FASTQ records with gzip compression. The writer is
multi-threaded and the actual write function call happens in another thread in
order to hide the IO latency.

Higher-level API
----------------
Apart from `KStream` class, this library provides another level of abstraction
which hides most details and provides very simple API on top of `KStream` for
working with sequence files: `SeqStreamIn` and `SeqStreamOut` for reading
and writing a sequence file respectively.  In order to prevent imposing any
unwanted external libraries (e.g. `zlib`) , the `SeqStream` class set are
defined in a separated header file (`seqio.h`) from the core library.

Reading a sequence file
-----------------------
These examples read FASTQ/A records one by one from either compressed or
uncompressed file.

Using `SeqStreamIn`:

```c++
#include <iostream>
#include "seqio.h"

using namespace klibpp;

int main(int argc, char* argv[])
{
  KSeq record;
  SeqStreamIn iss("file.dat");
  while (iss >> record) {
    std::cout << record.name << std::endl;
    if (!record.comment.empty()) std::cout << record.comment << std::endl;
    std::cout << record.seq << std::endl;
    if (!record.qual.empty()) std::cout << record.qual << std::endl;
  }
}
```

Using `KStream`:

```c++
#include <iostream>
#include <zlib>
#include "kseq++.h"

using namespace klibpp;

int main(int argc, char* argv[])
{
  KSeq record;
  gzFile fp = gzopen(filename, "r");
  auto ks = make_kstream(fp, gzread, mode::in);
  // auto ks = KStream(fp, gzread, mode::in);  // C++17
  // auto ks = KStreamIn(fp, gzread);  // C++17
  while (ks >> record) {
    std::cout << record.name << std::endl;
    if (!record.comment.empty()) std::cout << record.comment << std::endl;
    std::cout << record.seq << std::endl;
    if (!record.qual.empty()) std::cout << record.qual << std::endl;
  }
  gzclose(fp);
}
```

Or records can be fetched and stored in a `std::vector< KSeq >` in chunks.

Using `SeqStreamIn`:

```c++
#include <iostream>
#include "seqio.h"

using namespace klibpp;

int main(int argc, char* argv[])
{
  SeqStreamIn iss("file.dat");
  auto records = iss.read();
  // auto records = iss.read(100);  // read a chunk of 100 records
}
```

Using `KStream`:

```c++
#include <iostream>
#include <zlib>
#include "kseq++.h"

using namespace klibpp;

int main(int argc, char* argv[])
{
  gzFile fp = gzopen(filename, "r");
  auto ks = make_ikstream(fp, gzread);
  auto records = ks.read();  // fetch all the records
  // auto records = ks.read(100);  // read a chunk of 100 records
  gzclose(fp);
}
```

Writing a sequence file
-----------------------
These examples write FASTA/Q records to an uncompressed file.

Using `SeqStreamIn`:

```c++
#include <iostream>
#include "seqio.h"

using namespace klibpp;

int main(int argc, char* argv[])
{
  SeqStreamOut oss("file.dat");
  for (KSeq const& r : records) oss << r;
}
```

Using `KStream`:

```c++
#include <iostream>
#include <zlib>
#include "kseq++.h"

using namespace klibpp;

int main(int argc, char* argv[])
{
  int fd = open(filename, O_WRONLY);
  auto ks = make_kstream(fd, write, mode::out);
  // auto ks = KStreamOut(fd, write);  // C++ 17
  // ...
  for (KSeq const& r : records) ks << r;
  ks << kend;
  close(fd);
}
```

While writing a record to a file, sequence and quality scores can be wrapped at
a certain length. The default wrapping length is 60 bps and can be customised by
`KStream::set_wraplen` method.

---
**NOTE**

The buffer will be flushed to the file when the `KStream` object goes out of the
scope. Otherwise, `ks << kend` is required to be called before closing the file
to make sure that there is no data loss.

There is no need to write `kend` to the stream if using `SeqStreamOut`.

---

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

- CPU: Intel&reg; Xeon&reg; CPU E3-1241 v3 @ 3.50GHz, 4 cores, 8 threads
- RAM: DDR3 1600 MHz, 16352 MB
- HDD: Seagate Desktop HDD 500GB, 16MB Cache, SATA-3
- OS: Debian GNU/Linux 9.4 (stretch), Linux 4.9.91-1-amd64-smp
- Compiler: GCC 6.3.0, compiled with optimisation level 3 (`-O3`)

### Result

#### Reading all records

| file         |     kseq++ |   kseq |  SeqAn | kseq++/read\* | SeqAn/readRecords\*\* |
| :----------- | ---------: | -----: | -----: | ------------: | --------------------: |
| dataset_A.fa | **2.35 s** |  2.5 s | 2.92 s |    **3.52 s** |                4.94 s |
| dataset_B.fa | **2.66 s** |  2.8 s | 3.34 s |    **3.74 s** |                9.82 s |
| dataset_C.fq | **2.56 s** | 2.46 s | 2.66 s |    **4.56 s** |                11.8 s |

\* storing all records in `std::vector`.

\*\* storing all records in `seqan::StringSet< seqan::CharString >`.

#### Writing all records

| file         | kseq++/plain | kseq++/gzipped |  SeqAn/plain |
| :----------- | -----------: | -------------: | -----------: |
| dataset_A.fa |    **2.3 s** |      **866 s** |       2.29 s |
| dataset_B.fa |   **2.19 s** |      **849 s** |       2.33 s |
| dataset_C.fq |   **1.94 s** |      **365 s** |       2.24 s |
