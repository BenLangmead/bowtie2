 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: file_base.h,v 1.2 2009/02/19 01:51:23 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FILE_BASE_H
#define SEQAN_HEADER_FILE_BASE_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

	// To override the system's default temporary directory use the following:
	//#define SEQAN_DEFAULT_TMPDIR "/var/tmp"

	// To use direct I/O access define SEQAN_DIRECTIO (not completely tested yet)
	//#define SEQAN_DIRECTIO


/**
.Spec.Sync:
..cat:Files
..general:Class.File
..summary:File structure supporting synchronous input/output access.
..signature:File<Sync<> >
..remarks:This class suports pseudo-asynchronous access methods, i.e. the methods to initiate a I/O request return after request completion.
*/

	template <typename TSpec = void>
    struct Sync;

/**
.Spec.Async:
..cat:Files
..general:Class.File
..summary:File structure supporting synchronous and asynchronous input/output access.
..signature:File<Async<> >
*/

	template <typename TSpec = void>
    struct Async;


/**
.Class.File:
..cat:Input/Output
..summary:Represents a file.
..signature:File<TSpec>
..param.TSpec:The specializing type.
...default:$Async<>$, see @Spec.Async@.
*/

	template <typename TSpec = Async<> >
    class File;

/**
.Spec.Chained:
..cat:Files
..general:Class.File
..summary:Splits a large file into a chain of smaller files.
..signature:File<Chained<FileSize, TFile> >
..param.FileSize:The maximal split file size in byte.
...default:2^31-1 (~2GB)
..param.TFile:Underlying @Class.File@ type.
...default:$File<>$, see @Class.File@.
..remarks:This file type uses a chain of $TFile$ files, whose file sizes are at most $FileSize$ bytes.
Chained Files should be used for file systems or $TFile$ types that don't support large files (e.g. FAT32, C-style FILE*).
..remarks:The chain can be used as if it were one contiguous file.
*/

	// chained file's default filesize is 2gb-1byte (fat16 filesize limitation)
	template < __int64 _FileSize = ~(((__int64)1) << 63), typename TFile = File<> >
	struct Chained;

/**
.Spec.Striped:
..cat:Files
..general:Class.File
..summary:Stripes a file across multiple files.
..signature:File<Chained<FileCount, TFile> >
..param.FileCount:The number of files used for striping.
...default:2
..param.TFile:Underlying @Class.File@ type.
...default:$File<>$, see @Class.File@.
..remarks:This file type uses a software striping without redundance (see RAID0) to accelerate I/O access when using more than one disks.
..remarks:Striped files should only be used in @Class.Pool@s or external Strings as they only support block operations and no random accesses.
*/

	template < unsigned _FileCount = 2, typename TFile = File<> >
	struct Striped;

    enum FileOpenMode {
        OPEN_RDONLY     = 1,
        OPEN_WRONLY     = 2,
        OPEN_RDWR       = 3,
        OPEN_MASK       = 3,
        OPEN_CREATE     = 4,
        OPEN_APPEND     = 8,
        OPEN_ASYNC      = 16,
		OPEN_TEMPORARY	= 32,
		OPEN_QUIET		= 128
    };

	template <typename T>
	struct DefaultOpenMode {
		enum { VALUE = (OPEN_RDWR + OPEN_CREATE) | OPEN_APPEND };
	};

	template <typename T>
	struct DefaultOpenTempMode {
		enum { VALUE = OPEN_RDWR + OPEN_CREATE };
	};

    enum FileSeekMode {
        SEEK_BEGIN   = 0,
        SEEK_CURRENT = 1
#ifndef SEEK_END
      , SEEK_END     = 2
#endif
    };


    //////////////////////////////////////////////////////////////////////////////
    // result type of asynch. functions
    // you have to call release(aRequest<T>) after a finished *event based* transfer
	struct aDummyRequest {};

/**
.Class.aRequest:
..cat:Input/Output
..summary:Associated with an asynchronous I/O request.
..signature:aRequest<TFile>
..param.TFile:A File type.
..remarks:This structure is used to identify asynchronous requests after their initiation.
*/

    template < typename T >
    struct aRequest
    {
        typedef aDummyRequest Type;
    };
/*
    //////////////////////////////////////////////////////////////////////////////
    // event to represent asynchronous transfers
    // you can wait for it or test it
    template < typename T >
    struct aEvent
    {
        typedef DummyEvent Type;
    };

    ////////////////////////////////////////////////////////////////////////////////
    // callback hint parameter type
    // hint lets you recognize the finished asynch. transfer in your own callback routine
    template < typename T >
    struct aHint
    {
        typedef void Type;
    };

    //////////////////////////////////////////////////////////////////////////////
    // callback function interface
    template < typename T >
    struct aCallback
    {
        typedef void Type(aHint<T> *);
    };

    //////////////////////////////////////////////////////////////////////////////
    // file queue interface
    template < typename T >
    struct aQueue
    {
        typedef Nothing Type;
    };
*/

    //////////////////////////////////////////////////////////////////////////////
    // generic open/close interface

/**
.Function.open:
..summary:Opens a file.
..cat:Input/Output
..signature:open(file, fileName[, openMode])
..param.file:A File object.
...type:Class.File
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
*/

    template < typename TSpec >
    inline bool open(File<TSpec> &me, const char *fileName, int openMode) {
        return me.open(fileName, openMode);
    }

    template < typename TSpec >
    inline bool open(File<TSpec> &me, const char *fileName) {
		return open(me, fileName, DefaultOpenMode<File<TSpec> >::VALUE);
    }

/**
.Function.openTemp:
..summary:Opens a temporary file.
..cat:Input/Output
..signature:openTemp(file)
..param.file:A File object.
...type:Class.File
..remarks:After closing this file will automatically be deleted.
..remarks:The openmode (see @Function.open@) is $OPEN_RDWR | OPEN_CREATE$.
..returns:A $bool$ which is $true$ on success.
*/

    template < typename TSpec >
    inline bool openTemp(File<TSpec> &me) {
        return me.openTemp();
    }

    template < typename TSpec >
    inline bool openTemp(File<TSpec> &me, int openMode) {
        return me.openTemp(openMode);
    }

    template < typename File >
    inline void reopen(File &, int) {
	}

/**
.Function.close:
..cat:Input/Output
..summary:Closes a file.
..signature:close(file)
..param.file:A File object.
...type:Class.File
..returns:A $bool$ which is $true$ on success.
*/

    template < typename TSpec >
    inline bool close(File<TSpec> & me) {
        return me.close();
    }

    template < typename TSpec >
    inline unsigned sectorSize(File<TSpec> const & /*me*/) {
        return 4096;
    }


    //////////////////////////////////////////////////////////////////////////////
    // generic read(At)/write(At) interface

/**
.Function.read:
..cat:Input/Output
..summary:Loads records from a file.
..signature:read(file, memPtr, count)
..param.file:A File object.
...type:Class.File
..param.memPtr:A pointer to the first destination record in memory.
..param.count:The amount of records to be read.
..returns:A $bool$ which is $true$ on success.
..remarks:The records are read from the position pointed by the current file pointer (see @Function.seek@).
*/

	template < typename TSpec, typename TValue, typename TSize >
    inline bool read(File<TSpec> & me, TValue *memPtr, TSize const count) {
		return me.read(memPtr, count * sizeof(TValue));
    }

/**
.Function.write:
..cat:Input/Output
..summary:Saves records to a file.
..signature:write(file, memPtr, count)
..param.file:A File object.
...type:Class.File
..param.memPtr:A pointer to the first source record in memory.
..param.count:The amount of records to be written.
..returns:A $bool$ which is $true$ on success.
..remarks:The records are written at the position pointed by the current file pointer (see @Function.seek@).
*/

	template < typename TSpec, typename TValue, typename TSize >
    inline bool write(File<TSpec> & me, TValue const *memPtr, TSize const count) {
		return me.write(memPtr, count * sizeof(TValue));
    }

/**
.Function.readAt:
..summary:Loads records from a specific position in a file.
..cat:Input/Output
..signature:readAt(file, memPtr, count, fileOfs)
..param.file:A File object.
...type:Class.File
..param.memPtr:A pointer to the first destination record in memory.
..param.count:The amount of records to be read.
..param.fileOfs:The absolute file position in bytes measured from the beginning.
..returns:A $bool$ which is $true$ on success.
*/

    template < typename TFile, typename TValue, typename TSize, typename TPos >
    inline bool readAt(TFile & me, TValue *memPtr, TSize const count, TPos const fileOfs) {
		typedef typename Position<TFile>::Type pos_t;
		seek(me, (pos_t)fileOfs * (pos_t)sizeof(TValue));
		return read(me, memPtr, count);
    }

/**
.Function.writeAt:
..summary:Saves records to a specific position in a file.
..cat:Input/Output
..signature:writeAt(file, memPtr, count, fileOfs)
..param.file:A File object.
...type:Class.File
..param.memPtr:A pointer to the first source record in memory.
..param.count:The amount of records to be written.
..param.fileOfs:The absolute file position in bytes measured from the beginning.
..returns:A $bool$ which is $true$ on success.
*/

    template < typename TFile, typename TValue, typename TSize, typename TPos >
    inline bool writeAt(TFile & me, TValue const *memPtr, TSize const count, TPos const fileOfs) {
		typedef typename Position<TFile>::Type pos_t;
		seek(me, (pos_t)fileOfs * (pos_t)sizeof(TValue));
		return write(me, memPtr, count);
    }



    //////////////////////////////////////////////////////////////////////////////
    // generic seek/tell/size/resize interface

/**
.Function.seek:
..summary:Changes the current file pointer.
..cat:Input/Output
..signature:seek(file, fileOfs[, origin])
..param.file:A File object.
...type:Class.File
..param.fileOfs:A file offset measured in bytes relative to $origin$.
..param.origin:Selects the origin from where to calculate the new position.
...default:$SEEK_BEGIN$
...remarks:For $SEEK_BEGIN$, $SEEK_CURRENT$, or $SEEK_END$ the origin is the beginning, the current pointer, or the end of the file.
..returns:The new file position measured in bytes from the beginning.
*/

	template < typename TSpec, typename TPos >
    inline typename Position< File<TSpec> >::Type seek(File<TSpec> &me, TPos const fileOfs, int origin) {
		typedef typename Position< File<TSpec> >::Type TFilePos;
		TFilePos newOfs = me.seek(fileOfs, origin);
        #ifdef SEQAN_DEBUG_OR_TEST_
			if (origin == SEEK_BEGIN && newOfs != (TFilePos)fileOfs) {
				::std::cerr << "seek returned " << ::std::hex << newOfs << " instead of " << fileOfs << ::std::dec << ::std::endl;
			}
        #endif
        return newOfs;
    }

	template < typename TSpec, typename TPos >
    inline typename Position< File<TSpec> >::Type seek(File<TSpec> &me, TPos const fileOfs) {
		return seek(me, fileOfs, SEEK_BEGIN);
	}
/**
.Function.tell:
..summary:Gets the current file pointer.
..cat:Input/Output
..signature:tell(file)
..param.file:A File object.
...type:Class.File
..returns:The current file position measured in bytes from the beginning.
*/

    template < typename TSpec >
    inline typename Position< File<TSpec> >::Type tell(File<TSpec> &me) {
        return me.tell();
    }

/**
.Function.rewind:
..summary:Sets the current file pointer to the beginning.
..cat:Input/Output
..signature:rewind(file)
..param.file:A File object.
...type:Class.File
..remarks:Calls @Function.seek@$(file, 0)$ by default.
*/

    template < typename File >
    inline void rewind(File &me) {
		seek(me, 0);
    }

/**
.Function.size:
..summary:Gets the file size.
..cat:Input/Output
..signature:size(file)
..param.file:A File object.
...type:Class.File
..returns:The file size measured in bytes.
*/

    template < typename TSpec >
    inline typename Size<File<TSpec> >::Type size(File<TSpec> &me) {
        typename Size<File<TSpec> >::Type old_pos = tell(me);
        typename Size<File<TSpec> >::Type result = seek(me, 0, SEEK_END);
        seek(me, old_pos, SEEK_BEGIN);
        return result;
    }

/**
.Function.resize:
..cat:Input/Output
..signature:resize(file, new_length)
..param.file:A File object.
...type:Class.File
..param.new_length:The new file size measured in bytes.
*/

    template < typename TSpec, typename TSize >
    inline void resize(File<TSpec> &me, TSize new_length) {
        typename Size<File<TSpec> >::Type old_pos = tell(me);
        seek(me, new_length, SEEK_BEGIN);
        setEOF(me);
        seek(me, old_pos, SEEK_BEGIN);
    }

/**
.Function.setEOF:
..summary:Sets the file end to the current pointer.
..cat:Input/Output
..signature:setEOF(file)
..param.file:A File object.
...type:Class.File
*/

    template < typename TSpec >
    inline bool setEOF(File<TSpec> &/*me*/) {
		return true;
	}


    //////////////////////////////////////////////////////////////////////
    // Pseudo asynchronous Methods
    //////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////
    // callback based read/write
/*
    template < typename File, typename TValue, typename TSize,
               typename aCallback, typename aHint >
    inline typename aRequest<File>::Type
    aread(File & me, TValue *memPtr, TSize const count,
        aCallback* cb, aHint* hint)
    {
        result = read(me, memPtr, count);
        cb(hint);
        return NULL;
    }

    template < typename File, typename TValue, typename TSize,
               typename aCallback, typename aHint >
    inline typename aRequest<File>::Type
    awrite(File & me, TValue const *memPtr, TSize const count,
        aCallback* cb, aHint* hint)
    {
        write(me, memPtr, count);
        cb(hint);
        return NULL;
    }

    template < typename File, typename TValue, typename TSize, typename TPos,
               typename aCallback, typename aHint >
    inline typename aRequest<File>::Type
    areadAt(File & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aCallback* cb, aHint* hint)
    {
        readAt(me, memPtr, count, fileOfs);
        cb(hint);
        return NULL;
    }

    template < typename File, typename TValue, typename TSize, typename TPos,
               typename aCallback, typename aHint >
    inline typename aRequest<File>::Type
    awriteAt(File & me, TValue const *memPtr, TSize const count, TPos const fileOfs,
        aCallback* cb, aHint* hint)
    {
        result = writeAt(me, memPtr, count, fileOfs);
        cb(hint);
        return NULL;
    }


    //////////////////////////////////////////////////////////////////////
    // event based read/write

    template < typename File, typename TValue, typename TSize,
               typename aEvent >
    inline typename aRequest<File>::Type
    aread(File & me, TValue *memPtr, TSize const count,
        aEvent &event)
    {
        read(me, memPtr, count);
        event.signal();
        return NULL;
    }

    template < typename File, typename TValue, typename TSize,
               typename aEvent >
    inline typename aRequest<File>::Type
    awrite(File & me, TValue const *memPtr, TSize const count,
        aEvent &event)
    {
        write(me, memPtr, count);
        event.signal();
        return NULL;
    }

    template < typename File, typename TValue, typename TSize, typename TPos,
               typename aEvent >
    inline typename aRequest<File>::Type
    areadAt(File & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aEvent &event)
    {
        readAt(me, memPtr, count, fileOfs);
        event.signal();
        return NULL;
    }

    template < typename File, typename TValue, typename TSize, typename TPos,
               typename aEvent >
    inline typename aRequest<File>::Type
    awriteAt(File & me, TValue const *memPtr, TSize const count, TPos const fileOfs,
        aEvent &event)
    {
        writeAt(me, memPtr, count, fileOfs);
        event.signal();
        return NULL;
    }
*/

    //////////////////////////////////////////////////////////////////////
    // queue-less request based pseudo asychronous read/write

/**
.Function.areadAt:
..summary:Asynchronously loads records from a specific position in a file.
..cat:Input/Output
..signature:areadAt(file, memPtr, count, fileOfs, request)
..param.file:A File object.
...type:Class.File
..param.memPtr:A pointer to the first destination record in memory.
..param.count:The amount of records to be read.
..param.fileOfs:The absolute file position in bytes measured from the beginning.
..param.request:Reference to a structure that will be associated with this asynchronous request.
...type:Class.aRequest
..returns:A $bool$ which is $true$ on success.
*/

    template < typename File, typename TValue, typename TSize, typename TPos,
               typename aRequest >
    inline bool
	areadAt(File & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aRequest &request)
    {
        return readAt(me, memPtr, count, fileOfs);
    }

/**
.Function.awriteAt:
..summary:Asynchronously saves records to a specific position in a file.
..cat:Input/Output
..signature:awriteAt(file, memPtr, count, fileOfs, request)
..param.file:A File object.
...type:Class.File
..param.memPtr:A pointer to the first source record in memory.
..param.count:The amount of records to be written.
..param.fileOfs:The absolute file position in bytes measured from the beginning.
..param.request:Reference to a structure that will be associated with this asynchronous request.
...type:Class.aRequest
..returns:A $bool$ which is $true$ on success.
*/

    template < typename File, typename TValue, typename TSize, typename TPos,
               typename aRequest >
    inline bool
	awriteAt(File & me, TValue const *memPtr, TSize const count, TPos const fileOfs,
        aRequest &request)
    {
        return writeAt(me, memPtr, count, fileOfs);
    }


	//////////////////////////////////////////////////////////////////////
    // pseudo queue specific functions

/**
.Function.flush:
..summary:Waits for all open requests to complete.
..cat:Input/Output
..signature:flush(file)
..param.file:A File object.
...type:Class.File
..remarks:$flush$ returns after all pending requests are completed.
*/

    template < typename TSpec >
    inline void flush(File<TSpec> &) {
	}

/**
.Function.waitFor:
..summary:Waits for an asynchronous request to complete.
..cat:Input/Output
..signature:waitFor(request[, timeout_millis])
..param.request:Reference to an aRequest object.
...type:Class.aRequest
..param.timeout_millis:Timout value in milliseconds.
...remarks:A value of 0 can be used to test for completion without waiting.
...default:Infinity.
..returns:A $bool$ which is $true$ on completion and $false$ on timeout.
..remarks:$waitFor$ suspends the calling process until $request$ is completed or after $timeout_millis$ milliseconds.
*/

    inline bool waitFor(aDummyRequest &) {
		return true;
	}

	template < typename TTime >
    inline bool waitFor(aDummyRequest &, TTime) {
		return true;
	}

	// deprecated
	template < typename TSpec, typename aRequest >
    inline void release(File<TSpec> &, aRequest &) {
	}

/**
.Function.cancel:
..summary:Cancels an asynchronous request.
..cat:Input/Output
..signature:cancel(file, request)
..param.file:A File object.
...type:Class.File
..param.request:Reference to an aRequest object.
...type:Class.aRequest
..returns:A $bool$ which is $true$ on success.
*/

    template < typename TSpec, typename aRequest >
    inline bool cancel(File<TSpec> &, aRequest &) {
		return true;
	}


	// little helpers

	template <typename T1, typename T2> inline
	T1 enclosingBlocks(T1 _size, T2 _blockSize) {
		return (_size + _blockSize - 1) / _blockSize;
	}

	template <typename T1, typename T2> inline
	T1 alignSize(T1 _size, T2 _aligning) {
        if (_size < _aligning)
            return _aligning;
        else
		    return (_size / _aligning) * (T1)_aligning;
	}

}

#endif
