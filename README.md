kseq++
======
kseq++ is a C++11 re-implementation of [kseq](https://github.com/lh3/seqtk/blob/master/kseq.h)[.h](https://github.com/attractivechaos/klib/blob/master/kseq.h)
by [Heng Li](https://github.com/lh3). The goal for re-implementation of `kseq` is
providing better API and resource management while preserving its flexibility
and performance. Like original kseq, this parser is based on generic stream
buffer and works with different file types. However, instead of using C macros,
it uses C++ templates.

It inherits all features from kseq (quoting from kseq homepage):

> - Parse both FASTA and FASTQ format, and even a mixture of FASTA and FASTQ records in one file.
> - Seamlessly adapt to gzipped compressed file when used with zlib.
> - Support multi-line FASTQ.
> - Work on a stream with an internal stream buffer.

while additionally provides:
- simpler and more readable API
- RAII-style memory management

The library also comes with a **FASTA/Q writer**. Like reading, it can write
mixed multi-line FASTA and FASTQ records with _gzip compression_. The writer is
multi-threaded and the actual write function call happens in another thread in
order to hide the IO latency.

The RAII-style class `KStream` is the core class which handles input and output
streams. Each FASTA or FASTQ record will be stored in a `KSeq` object.

This library provides another layer of abstraction which hides most details and
provides very simple API on top of `KStream`: `SeqStreamIn` and `SeqStreamOut`
classes for reading and writing a sequence file respectively with exactly the
same interface. It is **highly recommended** to use these classes unless you
intent to use low-level interface like changing buffer size or use custom stream
type.

Looking for a quick start guide?
--------------------------------
Jump to [Examples](#examples) or [Installation](#installation).

KStream (`kseq++.hpp`)
----------------------
`KStream` is a generic, template class with the following template parameters
which are usually inferred by the compiler when constructed (so, there is no
need to provide them manually):
- `TFile`: type of the underlying stream/file (e.g. `gzFile`)
- `TFunc`: type of the read/write function corresponding to `TFile` (e.g.
  `int (*)(gzFile_s*, const void*, unsigned int)` for an output stream with
  `gzFile` as underlying file type)
- `TSpec`: stream opening mode (with values: `mode::in` or `mode::out`)

The template parameters are inferred by compiler in C++17 when instantiated by
calling their constructors. `make_kstream` function family also construct
`KStream`s which might be useful for inferring template parameters when using
older standards; e.g. C++11 or C++14.

To construct an instance, it requires at least three arguments: 1) the file
object/pointer/descriptor (can be of any type), 2) its corresponding read/write
function, and 3) stream opening mode (see [Examples](#examples)).

Higher-level API (`seqio.hpp`)
------------------------------
This header file defines `SeqStream` class set: i.e. `SeqStreamIn` and
`SeqStreamOut`. `SeqStream` classes are inherited from `KStream` with simpler
constructors using sensible defaults. They do not define any new method or
override inherited ones. So, they can be treated the same way as `KStream`.

In order to prevent imposing any unwanted external libraries (e.g. `zlib`) , the
`SeqStream` class set are defined in a separated header file (`seqio.hpp`) from
the core library.

Examples
--------

### Reading a sequence file
These examples read FASTQ/A records one by one from either compressed or
uncompressed file.

Using `SeqStreamIn`:

```c++
#include <iostream>
#include <kseq++/seqio.hpp>

using namespace klibpp;

int main(int argc, char* argv[])
{
  KSeq record;
  SeqStreamIn iss("file.fq.gz");
  while (iss >> record) {
    std::cout << record.name << std::endl;
    if (!record.comment.empty()) std::cout << record.comment << std::endl;
    std::cout << record.seq << std::endl;
    if (!record.qual.empty()) std::cout << record.qual << std::endl;
  }
}
```

<details>
<summary>Low-level API</summary>

Using `KStream`

```c++
#include <iostream>
#include <zlib>
#include <kseq++/kseq++.hpp>

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
</details>

Or records can be fetched and stored in a `std::vector< KSeq >` in chunks.

Using `SeqStreamIn`:

```c++
#include <iostream>
#include <kseq++/seqio.hpp>

using namespace klibpp;

int main(int argc, char* argv[])
{
  SeqStreamIn iss("file.fq");
  auto records = iss.read();
  // auto records = iss.read(100);  // read a chunk of 100 records
}
```

<details>
<summary>Low-level API</summary>

Using `KStream`

```c++
#include <iostream>
#include <zlib>
#include <kseq++/kseq++.hpp>

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
</details>

### Writing a sequence file
These examples write FASTA/Q records to an uncompressed file.

Using `SeqStreamIn`:

```c++
#include <iostream>
#include <kseq++/seqio.hpp>

using namespace klibpp;

int main(int argc, char* argv[])
{
  SeqStreamOut oss("file.dat");
  for (KSeq const& r : records) oss << r;
}
```

<details>
<summary>Low-level API</summary>

Using `KStream`

```c++
#include <iostream>
#include <zlib>
#include <kseq++/kseq++.hpp>

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
</details>

Another example for writing a series of FASTQ records to a gzipped file in
_FASTA_ format:

```c++
#include <iostream>
#include <kseq++/seqio.hpp>

using namespace klibpp;

int main(int argc, char* argv[])
{
  /* let `record` be a list of FASTQ records */
  SeqStreamOut oss("file.fa.gz", /* compression */ true, format::fasta);
  for (KSeq const& r : records) oss << r;
}
```

* * *
**NOTE**

The buffer will be flushed to the file when the `KStream` object goes out of the
scope. Otherwise, `ks << kend` is required to be called before closing the file
to make sure that there is no data loss.

There is no need to write `kend` to the stream if using `SeqStreamOut`.

* * *

#### Wrapping seq/qual lines

While writing a record to a file, sequence and quality scores can be wrapped at
a certain length. The default wrapping length for FASTA format is 60 bps and can
be customised by `KStream::set_wraplen` method. For FASTQ format -- i.e. when
the format is explicitly set to `format::fastq` -- output sequence and quality
string are not wrapped by default.

Wrapping can be disabled or enable by `KStream::set_nowrapping` and
`KStream::set_wrapping` methods respectively. The latter reset the wrapping
length to the default value (60 bps).

#### Formatting

The default behaviour is to write a record in FASTQ format if it has quality
information. Otherwise, i.e. when the quality string is empty, the record will
be written in FASTA format. So, the output might be a mixture of FASTQ and FASTA
records. However, the output format can be forced by using `format::fasta` and
`format::fastq` modifiers. For example:

```c++
out << format::fasta << fastq_record;
out << another_record;  // all other calls after this will also be in FASTA format.
```

will write a FASTQ record in FASTA format. These modifiers affect all writes
after them until another modifier is used. The `format::mix` modifier reverts
the behaviour to default.

* * *
**NOTE**

Writing a FASTA record in FASTQ format throws an exception unless the record is
empty (a record with empty sequence and quality string).

* * *

Installation
------------
kseq++ is a header-only library and can be simply included in a project. Use
the package provided in the
[Releases](https://github.com/cartoonist/kseqpp/releases) section and copy
`include/kseq++` to your project tree.

The `kseq++.hpp` is the core header file and `seqio.hpp` is optional and only
needs to be included when using higher-level API (see
[above](#higher-level-api-seqio.hpp)). The latter requires `zlib` as dependency
which should be linked.

There are also other ways to install the library:

### From source
Installing from source requires CMake>= 3.10:

``` shell
git clone https://github.com/cartoonist/kseqpp
cd kseqpp
mkdir build && cd build
cmake .. # -DCMAKE_INSTALL_PREFIX=/path/to/custom/install/prefix (optional)
make install
```

### From conda
It is also distributed on bioconda:
``` shell
conda install -c bioconda kseqpp
```

Development
-----------

### CMake integration
After installing the library, you can import the library to your project using
`find_package`. It imports `kseq++::kseq++` target which can be passed to
`target_include_directories` and `target_link_libraries` calls. This is a sample
CMake file for building `myprogram` which uses the library:

``` cmake
cmake_minimum_required(VERSION 3.10)
project(myprogram VERSION 0.0.1 LANGUAGES CXX)

find_package(kseq++ REQUIRED)

set(SOURCES "src/main.cpp")
add_executable(myprogram ${SOURCES})
target_include_directories(myprogram
  PRIVATE kseq++::kseq++)
target_link_libraries(myprogram
  PRIVATE kseq++::kseq++)
```

CMake options:
- for building tests: `-DBUILD_TESTING=on`
- for building benchmark: `-DBUILD_BENCHMARKING=on`

Benchmark
---------
**NOTE**: The results below are based on older versions of kseq++ and `kseq.h`.
  - [ ] TODO Update benchmark

**NOTE**: It is fair to say that kseq++ comes with a very negligible overhead
and is _almost_ as fast as `kseq.h` (in 'read' mode) with an idiomatic C++ API
and more convinient resource management. The original `kseq.h` does not support
writing FASTA/Q files.

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

### Result (for kseq++ v0.1.4)

#### Reading all records

| file         |     kseq++ |   kseq |  SeqAn | kseq++/read\* | SeqAn/readRecords\*\* |
| :----------- | ---------: | -----: | -----: | ------------: | --------------------: |
| dataset_A.fa | **2.35 s** |  2.5 s | 2.92 s |    **3.52 s** |                4.94 s |
| dataset_B.fa | **2.66 s** |  2.8 s | 3.34 s |    **3.74 s** |                9.82 s |
| dataset_C.fq | **2.56 s** | 2.46 s | 2.66 s |    **4.56 s** |                11.8 s |

\* storing all records in `std::vector`.

\*\* storing all records in `seqan2::StringSet< seqan2::CharString >`.

#### Writing all records

| file         | kseq++/plain | kseq++/gzipped |  SeqAn/plain |
| :----------- | -----------: | -------------: | -----------: |
| dataset_A.fa |    **2.3 s** |      **866 s** |       2.29 s |
| dataset_B.fa |   **2.19 s** |      **849 s** |       2.33 s |
| dataset_C.fq |   **1.94 s** |      **365 s** |       2.24 s |
