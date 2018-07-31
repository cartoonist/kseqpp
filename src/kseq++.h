/**
 *    @file  kseq++.h
 *   @brief  C++ implementation of kseq library.
 *
 *  This is a header-only library re-implementing the original kseq library.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Sun Jul 15, 2018  19:15
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2018, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef  KSEQPP_H__
#define  KSEQPP_H__

#include <cassert>
#include <cctype>
#include <cstring>
#include <cstdlib>
#include <ios>

// versioning
#define KLIBPP_MAJOR 0
#define KLIBPP_MINOR 0
#define KLIBPP_REVISION 1

namespace klibpp {
  struct KSeq {  // kseq_t
    std::string name;
    std::string comment;
    std::string seq;
    std::string qual;
    inline void clear( ) {
      name.clear();
      comment.clear();
      seq.clear();
      qual.clear();
    }
  };

  enum KStreamMode {
    in,
    out,
  };

  template< typename TFile,
            typename TFunc,
            int TBufSize = 16384 >
    class KStream {  // kstream_t
      protected:
        /* Separators */
        constexpr static int SEP_SPACE = 0;  // isspace(): \t, \n, \v, \f, \r
        constexpr static int SEP_TAB = 1;    // isspace() && !' '
        constexpr static int SEP_LINE = 2;   // line separator: "\n" (Unix) or "\r\n" (Windows)
        constexpr static int SEP_MAX = 2;
        /* Consts */
        constexpr static unsigned int DEFAULT_WRAPLEN = 60;
        /* Data members */
        unsigned char buf[ TBufSize ];       /**< @brief character buffer */
        int begin;                           /**< @brief begin buffer index */
        int end;                             /**< @brief end buffer index or error flag if -1 */
        bool is_eof;                         /**< @brief eof flag */
        bool is_tqs;                         /**< @brief truncated quality string flag */
        bool is_ready;                       /**< @brief next record ready flag */
        bool last;                           /**< @brief last read was successful */
        unsigned int wraplen;                /**< @brief line wrap length */
        KStreamMode mode;                    /**< @brief stream mode */
        TFile f;                             /**< @brief file handler */
        TFunc func;                          /**< @brief read/write function */
      public:
        KStream( TFile f_, TFunc func_, KStreamMode m_=KStreamMode::in )  // ks_init
          : wraplen( DEFAULT_WRAPLEN ), mode( m_ ), f( f_ ), func( func_ )
        {
          this->rewind();
        }

        KStream( KStream const& ) = delete;
        KStream& operator=( KStream const& ) = delete;
        KStream( KStream&& ) = default;
        KStream& operator=( KStream&& ) = default;
        ~KStream( ) = default;

        /* Mutators */
          inline void
        set_wraplen( unsigned int len )
        {
          this->wraplen = len;
        }
        /* Methods */
          inline bool
        err( ) const  // ks_err
        {
          return this->end == -1;
        }

          inline bool
        eof( ) const  // ks_eof
        {
          return this->is_eof && this->begin >= this->end;
        }

          inline bool
        tqs( ) const
        {
          return this->is_tqs;
        }

          inline bool
        fail( ) const
        {
          return this->err() || this->tqs() || ( this->eof() && !this->last );
        }

          inline void
        rewind( )  // ks_rewind
        {
          this->begin = 0;
          this->end = 0;
          this->is_eof = false;
          this->is_tqs = false;
          this->is_ready = false;
          this->last = false;
        }

          inline KStream&
        operator>>( KSeq& rec )  // kseq_read
        {
          int c;
          this->last = false;
          if ( !this->is_ready ) {  // then jump to the next header line
            while ( this->getc( c ) && c != '>' && c != '@' );
            if ( this->fail() ) return *this;
            this->is_ready = true;
          }  // else: the first header char has been read in the previous call
          rec.clear();  // reset all members
          if ( !this->getuntil( KStream::SEP_SPACE, rec.name, &c ) ) return *this;
          if ( c != '\n' ) {  // read FASTA/Q comment
            this->getuntil( KStream::SEP_LINE, rec.comment, nullptr );
          }
          while ( this->getc( c ) && c != '>' && c != '@' && c != '+' ) {
            if ( c == '\n' ) continue;  // skip empty lines
            rec.seq += c;
            this->getuntil( KStream::SEP_LINE, rec.seq, nullptr, true ); // read the rest of the line
          }
          this->last = true;
          if ( c == '>' || c == '@' ) this->is_ready = true;  // the first header char has been read
          if ( c != '+' ) return *this;  // FASTA
          while ( this->getc( c ) && c != '\n' );  // skip the rest of '+' line
          if ( this->eof() ) {  // error: no quality string
            this->is_tqs = true;
            return *this;
          }
          while ( this->getuntil( KStream::SEP_LINE, rec.qual, nullptr, true ) &&
              rec.qual.size() < rec.seq.size() );
          if ( this->err() ) return *this;
          this->is_ready = false;  // we have not come to the next header line
          if ( rec.seq.size() != rec.qual.size() ) {  // error: qual string is of a different length
            this->is_tqs = true;  // should return here
          }

          return *this;
        }

          inline KStream&
        operator<<( const KSeq& rec )
        {
          if ( rec.qual.empty() ) this->puts( ">" );  // FASTA record
          else this->puts( "@" );  // FASTQ record
          this->puts( rec.name );
          if ( !rec.comment.empty() ) {
            this->puts( " " );
            this->puts( rec.comment );
          }
          this->puts( "\n" );
          this->puts( rec.seq, true );
          if ( !rec.qual.empty() ) {
            this->puts( "\n+\n" );
            this->puts( rec.qual, true );
          }
          this->puts( "\n" );
          return *this;
        }

        operator bool( ) const
        {
          return !this->fail();
        }
        /* Low-level methods */
          inline bool
        getc( int& c ) noexcept  // ks_getc
        {
          // ready
          if ( this->begin < this->end ) {
            c = this->_nextc();
            return true;
          }
          // error
          if ( this->err() || this->eof() ) return false;
          // fetch
          this->begin = 0;
          this->end = this->func( this->f, this->buf, TBufSize );
          if ( this->end <= 0 ) {  // err if end == -1 and eof if 0
            this->is_eof = true;
            return false;
          }
          c = this->_nextc();
          return true;
        }

          inline bool
        puts( std::string const& s, bool wrap=false ) noexcept
        {
          if ( this->err() ) return false;

          std::string::size_type cursor = 0;
          std::string::size_type len = 0;
          while ( cursor != s.size() ) {
            assert( cursor < s.size() );
            len = s.size() - cursor;
            if ( wrap && this->wraplen < len ) len = this->wraplen;
            if ( wrap && cursor != 0 ) {
              if ( this->func( this->f, "\n", 1 ) == 0 ) {
                this->end = -1;
                break;
              }
            }
            if ( this->func( this->f, &s[ cursor ], len ) == 0 ) {
              this->end = -1;
              break;
            }
            cursor += len;
          }
          return !this->err();
        }

          inline bool
        getuntil( int delimiter, std::string& str, int *dret, bool append=false )  // ks_getuntil
          noexcept
        {
          int c;
          bool gotany = false;
          if ( dret ) *dret = 0;
          if ( !append ) str.clear();
          int i = -1;
          do {
            if ( !this->getc( c ) ) break;
            --this->begin;
            if ( delimiter == KStream::SEP_LINE ) {
              for ( i = this->begin; i < this->end; ++i ) {
                if ( this->buf[ i ] == '\n' ) break;
              }
            }
            else if ( delimiter > KStream::SEP_MAX ) {
              for ( i = this->begin; i < this->end; ++i ) {
                if ( this->buf[ i ] == delimiter ) break;
              }
            }
            else if ( delimiter == KStream::SEP_SPACE ) {
              for ( i = this->begin; i < this->end; ++i ) {
                if ( std::isspace( this->buf[ i ] ) ) break;
              }
            }
            else if ( delimiter == KStream::SEP_TAB ) {
              for ( i = this->begin; i < this->end; ++i ) {
                if ( std::isspace( this->buf[ i ] ) && this->buf[ i ] != ' ' ) break;
              }
            }
            else {
              assert( false );  // it should not reach here
              return false;  // when assert is replaced by NOOP
            }

            gotany = true;
            str.resize( str.size() + i - this->begin );
            std::copy( this->buf + this->begin, this->buf + i, str.end() - i + this->begin );
            this->begin = i + 1;
          } while ( i >= this->end );

          if ( this->err() || ( this->eof() && !gotany ) ) return false;

          assert( i != -1 );
          if ( !this->eof() && dret ) *dret = this->buf[ i ];
          if ( delimiter == KStream::SEP_LINE && !str.empty() && str.back() == '\r' ) {
            str.pop_back();
          }
          return true;
        }
      private:
        /* Methods */
          inline int
        _nextc( ) noexcept
        {
          return static_cast< int >( this->buf[ this->begin++ ] );
        }
    };

  template< typename TFile, typename TFunc, typename... Args >
      inline KStream< std::decay_t< TFile >, std::decay_t< TFunc > >
    make_kstream( TFile&& file, TFunc&& func, Args&&... args )
    {
      return KStream< std::decay_t< TFile >, std::decay_t< TFunc > >(
          std::forward< TFile >( file ), std::forward< TFunc >( func ),
          std::forward< Args >( args )... );
    }

  template< int TBufSize, typename TFile, typename TFunc, typename... Args >
      inline KStream< std::decay_t< TFile >, std::decay_t< TFunc >, TBufSize >
    make_kstream( TFile&& file, TFunc&& func, Args&&... args )
    {
      return KStream< std::decay_t< TFile >, std::decay_t< TFunc >, TBufSize >(
          std::forward< TFile >( file ), std::forward< TFunc >( func ),
          std::forward< Args >( args )... );
    }
}  /* -----  end of namespace klibpp  ----- */
#endif  /* ----- #ifndef KSEQPP_H__  ----- */
