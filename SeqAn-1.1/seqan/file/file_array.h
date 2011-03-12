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
  $Id: file_array.h,v 1.1 2008/08/25 16:20:03 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FILE_ARRAY_H
#define SEQAN_HEADER_FILE_ARRAY_H

#include <sstream>
#include <iomanip>


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

	//template < __int64 _FileSize = 2*1024*1024*1024-1, typename TFile = File<> >
	//struct Chained;

	//template < unsigned _FileCount = 2, typename TFile = File<> >
	//struct Striped;


    template < __int64 _FileSize, typename TFile >
    struct Size< File< Chained<_FileSize, TFile> > >
    {
        typedef __int64 Type;
    };

    template < __int64 _FileSize, typename TFile >
    struct Position< File< Chained<_FileSize, TFile> > >
    {
        typedef __int64 Type;
    };

    template < __int64 _FileSize, typename TFile >
    struct Difference< File< Chained<_FileSize, TFile> > >
    {
        typedef __int64 Type;
    };

    template < __int64 _FileSize, typename TFile >
    struct aRequest< File< Chained<_FileSize, TFile> > >
    {
		typedef typename aRequest<TFile>::Type Type;
    };


    template < unsigned _FileCount, typename TFile >
    struct Size< File< Striped<_FileCount, TFile> > >
    {
        typedef __int64 Type;
    };

    template < unsigned _FileCount, typename TFile >
    struct Position< File< Striped<_FileCount, TFile> > >
    {
        typedef __int64 Type;
    };

    template < unsigned _FileCount, typename TFile >
    struct Difference< File< Striped<_FileCount, TFile> > >
    {
        typedef __int64 Type;
    };

    template < unsigned _FileCount, typename TFile >
    struct aRequest< File< Striped<_FileCount, TFile> > >
    {
		typedef typename aRequest<TFile>::Type Type;
    };


	template < unsigned _FileCount, typename TFile >
	class File< Striped<_FileCount, TFile> >: public Tuple< TFile, _FileCount > {
		File(void *dummy = NULL) {}	// to be compatible with the FILE*(NULL) constructor
		operator bool() const { return (*this)[0]; }
	};

    template < __int64 _FileSize, typename TFile >
	class File< Chained<_FileSize, TFile> >: public String< TFile > {
		typedef String< TFile > Base;

		::std::string	baseName;
		int				openMode;
		__int64			fileSize;
		bool			temporary;

		File(void *dummy = NULL) :	// to be compatible with the FILE*(NULL) constructor
			fileSize(0),
			_realign(false) {}

	private:
		
		bool _realign;
	
		template < typename TSize, typename TValue >
		inline void _alignFloor(TSize _size, TValue const *) {
			__int64 alignment = sizeof(TValue) * sectorSize(TFile());
			fileSize = (_size / alignment) * alignment;
		}

		template < typename TSize, typename TValue >
		inline void _alignCeil(TSize _size, TValue const *) {
			__int64 alignment = sizeof(TValue) * sectorSize(TFile());
			fileSize = ((_size + alignment - 1) / alignment) * alignment;
		}

	public:
	
		inline ::std::string getFileName(int i) const { 
			::std::stringstream strm;
			strm << baseName << '.' << ::std::setfill('0') << ::std::setw(3) << i;
			return strm.str();
		}

		inline operator bool() const { 
			return (*this)[0]; 
		}

		inline unsigned fileCount() const {
			return length(*(Base*)this);
		}

		inline TFile& getFile(int fileNo) {
			unsigned _oldFileCount = fileCount();
			if (fileNo >= _oldFileCount) {
				resize(*(Base*)this, fileNo + 1);
				for(unsigned i = _oldFileCount; i <= fileNo; ++i)
					if (temporary)
						openTemp((*this)[i], openMode);
					else
						open((*this)[i], getFileName(i).c_str(), openMode);
			}
			return (*this)[fileNo];
		}

		inline void tryOpen() {
			unsigned fileCount = 0;
			while (fileExists(getFileName(fileCount).c_str())) ++fileCount;
			if (fileCount) {
				fileSize = size(getFile(0));
				_realign = (fileCount == 1);
				getFile(fileCount - 1);
			} 
		}

		// fileSize dependent functions

		template < typename TValue >
		inline void adjustFileSize(TValue const *dummy) {
			if (_realign) {
				_alignCeil(fileSize, dummy);
				_realign = false;
				if (fileSize < _FileSize)
					fileSize = 0;
			}
			if (!fileSize) 
				_alignFloor(_FileSize, dummy);
		}

		template < typename TPos, typename TOffset, typename TValue >
		inline TFile& getFileAndOffset(TPos offset, TOffset &fileOffset, TValue const *dummy) {
			adjustFileSize(dummy);
			offset *= sizeof(TValue);
			fileOffset = (offset % fileSize) / sizeof(TValue);
			return getFile(offset / fileSize);
		}

		template < typename TOffset, typename TValue >
		inline __int64 restAt(TOffset fileOffset, TValue const *dummy) {
			adjustFileSize(dummy);
			__int64 restBytes = fileSize;
			restBytes -= fileOffset * sizeof(TValue);
			return restBytes / sizeof(TValue);
		}

		inline void resizeArray(__int64 _newSize) {
			if (fileSize) {
				unsigned _oldFileCount = fileCount();
				unsigned _newFileCount = enclosingBlocks(_newSize, fileSize);
				for(unsigned i = _newFileCount; i < _oldFileCount; ++i) {
					close((*this)[i]);
					if (!temporary) fileUnlink(getFileName(i).c_str());
				}
				resize(*(Base*)this, _newFileCount);
				if (_newFileCount) {
					typename Size<TFile>::Type lastFileSize = _newSize % fileSize;
					if (fileSize) resize((*this)[_newFileCount - 1], lastFileSize);
				}
			}
		}

        inline void clearInternals() {
			clear(*(Base*)this);
            fileSize = 0;
            _realign = false;
		}
	};


    //////////////////////////////////////////////////////////////////////////////
    // generic open/close interface
    template < typename TFileArray >
    inline bool _openTempFArray(TFileArray &me, int openMode) {
		bool result = true;
		for(int i = 0; i < length(me); ++i)
			result &= openTemp(me[i], openMode);
		return result;
    }

    template < typename TFileArray >
    inline bool _openTempFArray(TFileArray &me) {
		return _openTempFArray(me, DefaultOpenTempMode<TFileArray>::VALUE);
	}

    template < typename TFileArray >
    inline bool _reopenFArray(TFileArray &me, int openMode) {
		bool result = true;
		for(int i = 0; i < length(me); ++i)
			result &= reopen(me[i], openMode);
		return result;
    }

    template < typename TFileArray >
    inline bool _closeFArray(TFileArray &me) {
		bool result = true;
		for(int i = 0; i < length(me); ++i)
			if (me[i]) result &= close(me[i]);
		return result;
    }

    template < typename TFileArray >
    inline unsigned _sectorSizeFArray(TFileArray &me, int openMode) {
		return sectorSize(me[0]);
    }

    template < typename TFileArray >
    inline typename Size<TFileArray>::Type
	_sizeFArray(TFileArray &me) {
        typename Size<TFileArray>::Type sum = 0;
		for(int i = 0; i < length(me); ++i)
			sum += size(me[i]);
		return sum;
    }

    template < typename TFileArray >
    inline bool _flushFArray(TFileArray &me) {
		bool result = true;
		for(int i = 0; i < length(me); ++i)
			result &= flush(me[i]);
		return result;
    }

    template < typename TFileArray, typename TRequest >
    inline bool _cancelFArray(TFileArray &me, TRequest &request) {
		bool result = true;
		for(int i = 0; i < length(me); ++i)
			result &= cancel(me[i], &request);
		return result;
    }


    //////////////////////////////////////////////////////////////////////////////
    // standard file array wrappers

    template < __int64 _FileSize, typename TFile >
	inline unsigned length(File< Chained<_FileSize, TFile> > const &me) {
		return me.fileCount();
	}

    template < unsigned _FileCount, typename TFile >
	inline unsigned length(File< Striped<_FileCount, TFile> > const &me) {
		return _FileCount;
	}

    template < __int64 _FileSize, typename TFile >
	inline bool open(File< Chained<_FileSize, TFile> > &me, const char *fileName, int openMode) {
		me.baseName = fileName;
		me.openMode = openMode;
		me.temporary = false;
		me.tryOpen();
		return true;
	}

    template < __int64 _FileSize, typename TFile >
	inline bool openTemp(File< Chained<_FileSize, TFile> > &me, int openMode) {
		me.openMode = openMode;
		me.temporary = true;
		return true;
	}

    template < unsigned _FileCount, typename TFile >
	inline bool openTemp(File< Striped<_FileCount, TFile> > &me, int openMode) {
		return _openTempFArray(me, openMode);
	}

    template < __int64 _FileSize, typename TFile >
	inline bool close(File< Chained<_FileSize, TFile> > &me) {
        _closeFArray(me);
        me.clearInternals();
        return true;
    }

    template < unsigned _FileCount, typename TFile >
	inline bool close(File< Striped<_FileCount, TFile> > &me) {	return _closeFArray(me); }

    template < __int64 _FileSize, typename TFile >
	__int64 size(File< Chained<_FileSize, TFile> > &me) {
		return _sizeFArray(me);
	}

    template < unsigned _FileCount, typename TFile >
	__int64 size(File< Striped<_FileCount, TFile> > &me) {
		return _sizeFArray(me);
	}

    template < __int64 _FileSize, typename TFile, typename TSize >
    inline void resize(File< Chained<_FileSize, TFile> > &me, TSize new_length) {
		me.resizeArray(new_length);
    }

    template < __int64 _FileSize, typename TFile, typename TValue, typename TSize >
	inline void allocate(File< Chained<_FileSize, TFile> > const &me, TValue* &data, TSize count) {
		allocate(me[0], data, count);
	}

    template < __int64 _FileSize, typename TFile, typename TValue, typename TSize >
	inline void deallocate(File< Chained<_FileSize, TFile> > const &me, TValue* &data, TSize count) {
		deallocate(me[0], data, count);
	}

    template < unsigned _FileCount, typename TFile, typename TValue, typename TSize >
	inline void allocate(File< Striped<_FileCount, TFile> > const &me, TValue* &data, TSize count) {
		allocate(me[0], data, count);
	}

    template < unsigned _FileCount, typename TFile, typename TValue, typename TSize >
	inline void deallocate(File< Striped<_FileCount, TFile> > const &me, TValue* &data, TSize count) {
		deallocate(me[0], data, count);
	}


    //////////////////////////////////////////////////////////////////////////////
    // read/write wrappers

    template < __int64 _FileSize, typename TFile, typename TValue, typename TSize, typename TOffset >
    inline bool readAt(File< Chained<_FileSize, TFile> > &me, TValue *memPtr, TSize count, TOffset offset) {
		TOffset fileOfs = 0;
		while (count) {
			TFile &file = me.getFileAndOffset(offset, fileOfs, memPtr);
			TSize xmitSize = _min(me.restAt(fileOfs, memPtr), (__int64)count);
			if (!readAt(file, memPtr, xmitSize, fileOfs)) return false;
			count -= xmitSize;
			offset += xmitSize;
			memPtr += xmitSize;
		}
		return true;
    }
    
    template < __int64 _FileSize, typename TFile, typename TValue, typename TSize, typename TOffset >
    inline bool writeAt(File< Chained<_FileSize, TFile> > &me, TValue const *memPtr, TSize count, TOffset offset) {
		TOffset fileOfs = 0;
		while (count) {
			TFile &file = me.getFileAndOffset(offset, fileOfs, memPtr);
			TSize xmitSize = _min(me.restAt(fileOfs, memPtr), (__int64)count);
			if (!writeAt(file, memPtr, xmitSize, fileOfs)) return false;
			count -= xmitSize;
			offset += xmitSize;
			memPtr += xmitSize;
		}
		return true;
    }

    template < __int64 _FileSize, typename TFile, typename TValue, typename TSize, typename TOffset, typename TRequest >
    inline bool areadAt(File< Chained<_FileSize, TFile> > &me, TValue *memPtr, TSize count, TOffset offset, TRequest &req) {
		TOffset fileOfs = 0;
		while (count) {
			TFile &file = me.getFileAndOffset(offset, fileOfs, memPtr);
			TSize xmitSize = _min(me.restAt(fileOfs, memPtr), (__int64)count);
			if (count != xmitSize) {
				if (!readAt(file, memPtr, xmitSize, fileOfs)) return false;
			} else
				if (!areadAt(file, memPtr, xmitSize, fileOfs, req)) return false;
			count -= xmitSize;
			offset += xmitSize;
			memPtr += xmitSize;
		}
		return true;
    }
    
    template < __int64 _FileSize, typename TFile, typename TValue, typename TSize, typename TOffset, typename TRequest  >
    inline bool awriteAt(File< Chained<_FileSize, TFile> > &me, TValue const *memPtr, TSize count, TOffset offset, TRequest &req) {
		TOffset fileOfs = 0;
		while (count) {
			TFile &file = me.getFileAndOffset(offset, fileOfs, memPtr);
			TSize xmitSize = _min(me.restAt(fileOfs, memPtr), (__int64)count);
			if (count != xmitSize) {
				if (!writeAt(file, memPtr, xmitSize, fileOfs)) return false;
			} else
				if (!awriteAt(file, memPtr, xmitSize, fileOfs, req)) return false;
			count -= xmitSize;
			offset += xmitSize;
			memPtr += xmitSize;
		}
		return true;
    }

}

#endif
