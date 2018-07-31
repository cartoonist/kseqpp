/**
 *    @file  kseq++_bench.cc
 *   @brief  Benchmark for kseq++.
 *
 *  Benchmark for kseq++ header file.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Tue Jul 17, 2018  23:02
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2018, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <zlib.h>
#include <iostream>
#include <ios>
#include <iomanip>
#include <ctime>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <fcntl.h>

#include "kseq++.h"
#include "kseq.h"
#include "seqan/seq_io.h"

KSEQ_INIT(gzFile, gzread)

#define SMALL_BUF_SIZE 4096
#define BIG_BUF_SIZE 65536


using namespace klibpp;

  int
main( int argc, char* argv[] )
{
  if ( argc == 1 ) {
    std::cerr << "Usage: " << argv[0] << " FILE" << std::endl;
    return EXIT_FAILURE;
  }

  gzFile fp;
  clock_t t;
  {
    std::string buf;
    buf.resize( SMALL_BUF_SIZE );
    fp = gzopen( argv[1], "r" );
    t = clock();
    while ( gzread( fp, &buf[0], SMALL_BUF_SIZE ) > 0 );
    auto d = static_cast< float >( clock() - t ) / CLOCKS_PER_SEC;
    std::cerr << "[gzread] " << std::setprecision(3) << d << " sec" << std::endl;
    gzclose( fp );
  }
  {
    fp = gzopen( argv[1], "r" );
    auto ks = make_kstream< SMALL_BUF_SIZE >( fp, gzread );
    KSeq record;
    int c;
    t = clock();
    while ( ks.getc( c ) );
    auto d = static_cast< float >( clock() - t ) / CLOCKS_PER_SEC;
    std::cerr << "[ks_getc] " << std::setprecision(3) << d << " sec" << std::endl;
    gzclose( fp );
  }
  {
    fp = gzopen( argv[1], "r" );
    auto ks = make_kstream< SMALL_BUF_SIZE >( fp, gzread );
    std::string s;
    int dret;
    t = clock();
    while ( ks.getuntil( '\n', s, &dret ) );
    auto d = static_cast< float >( clock() - t ) / CLOCKS_PER_SEC;
    std::cerr << "[ks_getuntil] " << std::setprecision(3) << d << " sec" << std::endl;
    gzclose( fp );
  }
  {
    fp = gzopen( argv[1], "r" );
    t = clock();
    while ( gzgetc( fp ) >= 0 );
    auto d = static_cast< float >( clock() - t ) / CLOCKS_PER_SEC;
    std::cerr << "[gzgetc] " << std::setprecision(3) << d << " sec" << std::endl;
    gzclose( fp );
  }
  {
    fp = gzopen( argv[1], "r" );
    std::string buf;
    buf.resize( SMALL_BUF_SIZE );
    t = clock();
    while ( gzgets( fp, &buf[0], SMALL_BUF_SIZE ) != Z_NULL );
    auto d = static_cast< float >( clock() - t ) / CLOCKS_PER_SEC;
    std::cerr << "[gzgets] " << std::setprecision(3) << d << " sec" << std::endl;
    gzclose( fp );
  }
  {
    FILE *fp;
    std::string s;
    t = clock();
    s.resize( BIG_BUF_SIZE );
    fp = fopen( argv[1], "r" );
    while ( fgets( &s[0], BIG_BUF_SIZE, fp ) );
    auto d = static_cast< float >( clock() - t ) / CLOCKS_PER_SEC;
    std::cerr << "[fgets] " << std::setprecision(3) << d << " sec" << std::endl;
    fclose( fp );
  }
  {
    int fd, dret;
    std::string s;
    t = clock();
    fd = open( argv[1], O_RDONLY );
    auto ks = make_kstream< BIG_BUF_SIZE >( fd, read );
    while ( ks.getuntil( '\n', s, &dret ) );
    auto d = static_cast< float >( clock() - t ) / CLOCKS_PER_SEC;
    std::cerr << "[kstream] " << std::setprecision(3) << d << " sec" << std::endl;
    close( fd );
  }
  {
    seqan::SeqFileIn f;
    seqan::CharString name;
    seqan::CharString str;
    seqan::CharString qual;
    t = clock();
    open( f, argv[1] );
    while ( !atEnd( f ) ) readRecord( name, str, qual, f );
    auto d = static_cast< float >( clock() - t ) / CLOCKS_PER_SEC;
    std::cerr << "[seqan] " << std::setprecision(3) << d << " sec" << std::endl;
    close( f );
  }
  {
    kseq_t *seq;
    int l;
    t = clock();
    fp = gzopen( argv[1], "r" );
    seq = kseq_init( fp );
    while ( ( l = kseq_read( seq ) ) >= 0 );
    auto d = static_cast< float >( clock() - t ) / CLOCKS_PER_SEC;
    std::cerr << "[kseq] " << std::setprecision(3) << d << " sec" << std::endl;
    kseq_destroy(seq);
    gzclose( fp );
  }
  {
    KSeq record;
    t = clock();
    fp = gzopen( argv[1], "r" );
    auto ks = make_kstream( fp, gzread );
    while (ks >> record);
    auto d = static_cast< float >( clock() - t ) / CLOCKS_PER_SEC;
    std::cerr << "[kseq++] " << std::setprecision(3) << d << " sec" << std::endl;
    gzclose( fp );
  }

  return EXIT_SUCCESS;
}
