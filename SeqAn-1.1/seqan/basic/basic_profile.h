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
  $Id: basic_profile.h,v 1.1 2008/08/25 16:20:02 langmead Exp $
 ==========================================================================*/

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_BASIC_PROFILE_H
#define SEQAN_HEADER_BASIC_PROFILE_H

// todo: substitute defines with inlines
#ifndef SEQAN_PROFILE

    #define SEQAN_PROSET(i,v)
    #define SEQAN_PROADD(i,v)
    #define SEQAN_PROSUB(i,v)
	#define SEQAN_PROVAL(i)				0
    #define SEQAN_PROEXTRAS(i)
    #define SEQAN_PROMARK(m)
    #define SEQAN_PROENDMARK(m)
    #define SEQAN_PRORESET
	#define SEQAN_PROGETTIME			0
    #define SEQAN_PROTIMESTART(a)
    #define SEQAN_PROTIMEDIFF(a)		0
	// replace malloc and free in external tools
	// with SEQAN_PROMALLOC and SEQAN_PROFREE to profile
	// their memory usage
    #define SEQAN_PROMALLOC(s) 			malloc(s)
    #define SEQAN_PROFREE(p) 			free(p)

#else

    #define SEQAN_PROSET(i,v)			_proSet(i,v)
    #define SEQAN_PROADD(i,v)			_proAdd(i,v)
    #define SEQAN_PROSUB(i,v)			_proSub(i,v)
	#define SEQAN_PROVAL(i)				(_proData<>::_proValue[i])
    #define SEQAN_PROEXTRAS(i)			{_proData<>::_proExtraCount = i;}
    #define SEQAN_PROMARK(m)			_proMark(m)
    #define SEQAN_PROENDMARK(m)			_proEndMark(m)
    #define SEQAN_PRORESET				_proReset()
	#define SEQAN_PROGETTIME			sysTime()
    #define SEQAN_PROTIMESTART(a)		_proFloat a = sysTime()
    #define SEQAN_PROTIMEDIFF(a)		(sysTime() - a)
    #define SEQAN_PROMALLOC(s)			_proMalloc(s)
    #define SEQAN_PROFREE(p)			_proFree(p)

#endif

#ifdef PLATFORM_WINDOWS
    typedef __int64   _proInt;
#else
    typedef int64_t _proInt;
#endif

    typedef double    _proFloat;


    typedef _proFloat _proTValue;

    enum _proConsts {
        SEQAN_PROPAGESIZE         = 4096, // B in byte
        SEQAN_PROFLOAT            = 0,
        SEQAN_PROINT              = 1,
        SEQAN_PROTIME             = 2,
        SEQAN_PROTYPEMASK         = 3,
        SEQAN_PROSTATE            = 4
    };

    enum _proValueIndex {
		SEQAN_PROSYSTIME		  = 0,
		SEQAN_PROCPUTIME		  = 1,
        SEQAN_PROMEMORY           = 2,    // current memory usage (state value)
        SEQAN_PROIO               = 3,    // IOs done (measured in Blocks of size B)
        SEQAN_PROIORANDOM         = 4,    // IOs calls done (read/write calls done)
        SEQAN_PROIOVOLUME         = 5,    // current disk usage (state value)
        SEQAN_PRODEPTH            = 6,    // algorithmic rec. depth or loop count
		SEQAN_PROOPENFILES		  = 7,	  // currently opened files
        SEQAN_PROIWAIT            = 8,    // waiting time (initiating)
        SEQAN_PROCWAIT            = 9,    // waiting time (completing)
		SEQAN_PROEXTRA1           = 10,
		SEQAN_PROEXTRA2           = 11,
		SEQAN_PROEXTRA3           = 12,
		SEQAN_PROINDEXCOUNT       = 13,
		SEQAN_PROEXTRACOUNT       = 3
    };

    const char _proValueType[] = {
		SEQAN_PROTIME, 
		SEQAN_PROTIME, 
        SEQAN_PROINT + SEQAN_PROSTATE, 
        SEQAN_PROINT,
        SEQAN_PROINT,
        SEQAN_PROINT + SEQAN_PROSTATE, 
        SEQAN_PROINT + SEQAN_PROSTATE, 
        SEQAN_PROINT + SEQAN_PROSTATE, 
        SEQAN_PROFLOAT,
        SEQAN_PROFLOAT,
        SEQAN_PROFLOAT + SEQAN_PROSTATE,
        SEQAN_PROFLOAT + SEQAN_PROSTATE,
        SEQAN_PROFLOAT + SEQAN_PROSTATE
    };

    typedef _proTValue _proTStates[SEQAN_PROINDEXCOUNT];
    typedef _proFloat  _proTTimes[SEQAN_PROINDEXCOUNT];



    struct _proFile;

	template <typename T = void>
	struct _proData
	{
		static _proTStates	_proValue;
		static _proTTimes	_proLastUpdate;
		static int			_proExtraCount;
	    
		static clock_t		_proCpuTimeLast;			// clock_t wraps around every 72mins
		static _proInt		_proCpuTimeOffset;			// we have to work around this

		static _proFile*	_proPFile;
		static _proFile*	_proPFileStream;
	};

	template <typename T> _proTStates	_proData<T>::_proValue = {};
	template <typename T> _proTStates	_proData<T>::_proLastUpdate = {};
	template <typename T> int			_proData<T>::_proExtraCount = 0;
	template <typename T> clock_t		_proData<T>::_proCpuTimeLast = 0;
	template <typename T> _proInt		_proData<T>::_proCpuTimeOffset = 0;
	template <typename T> _proFile*		_proData<T>::_proPFile = NULL;
	template <typename T> _proFile*		_proData<T>::_proPFileStream = NULL;


	inline _proFile* & _proPFile()			{ return _proData<>::_proPFile; }
	inline _proFile* & _proPFileStream()	{ return _proData<>::_proPFileStream; }


// HINT: The unit of all time functions is second.
    inline _proFloat cpuTime() {
    	clock_t now = clock();
    	if (_proData<>::_proCpuTimeLast > now) {		// test for time wrap
    		_proData<>::_proCpuTimeOffset += (~0u);		// got one
    		_proData<>::_proCpuTimeOffset ++;
//    		printf("\n!!WRAP!! old:%d, now:%d    ofs:%d\n",_proData<>::_proCpuTimeLast,now,_proData<>::_proCpuTimeOffset);
    	}
   		_proData<>::_proCpuTimeLast = now;
    	return (_proData<>::_proCpuTimeOffset + now) / (_proFloat)CLOCKS_PER_SEC;
   	}

    #ifdef PLATFORM_WINDOWS
//        inline _proFloat sysTime() { return GetTickCount() * 1e-3; }
		inline _proFloat sysTime() { return ( (_proFloat) clock() ) / CLOCKS_PER_SEC; }
    #else

		#include <unistd.h>
		#if _POSIX_TIMERS > 0
			#ifndef SEQAN_USE_CLOCKGETTIME
			#define SEQAN_USE_CLOCKGETTIME
			#endif
		#endif
		
		#ifndef SEQAN_USE_CLOCKGETTIME
		/* some systems e.g. darwin have no clock_gettime */
		
			#include <sys/time.h>
			
			inline _proFloat sysTime() {
				struct timeval tp;
				gettimeofday(&tp, NULL);
				return tp.tv_sec + tp.tv_usec * 1e-6;
			}

		#else

			inline _proFloat sysTime() {
			  /*struct timespec tp;
				clock_gettime(CLOCK_MONOTONIC, &tp);
				  return tp.tv_sec + tp.tv_nsec * 1e-9;*/
				return 0; // BTL: to compile under cygwin
			}

	    #endif

    #endif

    
    struct _proFile {

        FILE   *out;
        bool   running;

        _proFloat dumpStep;            // 0 .. manual dump mode, >0 .. live stream
        _proFloat dumpNext;        

        _proTStates all, last;
        ::std::string mark;
        unsigned	lines;

        _proFile() {
            running = false;
        }

        _proFile(char const *fname, _proFloat _dumpStep = 300.0) { // five minutes default dump interval
            running = false;
            start(fname, _dumpStep);
        }

        ~_proFile() {
            if (running) stop();
        }

        inline void start(char const *fname, _proFloat _dumpStep = 300.0, bool append = false) {
            if (append)
                out = fopen(fname, "a");
            else {
                out = fopen(fname, "w");
                dumpHeader();
            }

            if (!out) printf("WARNING: proFile could not be opened.\n");

			setTime(_proData<>::_proValue);
            syncAll(all);
            syncAll(last);
            running      = true;
            lines		 = 0;
            dumpStep     = _dumpStep;
            dumpNext     = sysTime();
            dump(last);
        }

        inline void stop() {
            dump(last);
            maximize(all, last);
            if (dumpStep == 0) {
                mark = "Zusammenfassung";
                dump(all);
            }
            fclose(out);
            running = false;
        }

        inline void syncTime(_proTStates &dst) {
            memcpy(dst, _proData<>::_proValue, 2 * sizeof(_proTValue));
        }

        inline void sync(_proTStates &dst) {
            memcpy(&(dst[2]), &(_proData<>::_proValue[2]), sizeof(_proTStates) - 2 * sizeof(_proTValue));
        }

        inline void syncAll(_proTStates &dst) {
            memcpy(dst, _proData<>::_proValue, sizeof(_proTStates));
        }

		inline static void setTime(_proTStates &dst) {
            dst[0] = sysTime();
			dst[1] = cpuTime();
		}

        inline void maximize(_proTStates &dst, _proTStates const &src) {
            for(int i = 0; i < SEQAN_PROINDEXCOUNT; ++i)
                if (((_proValueType[i] & SEQAN_PROSTATE) != 0))
                    if (dst[i] < src[i])
                        dst[i] = src[i];
        }

        inline void dumpTab() {
            if (!bol)
                fprintf(out, " \t");
            bol = false;
        }

        inline void dumpEndl() { fprintf(out, "\n"); }

        inline void dumpHeader() {
            fprintf(out, "\"Echtzeit\"\t\"CPU-Zeit\"\t\"Speicher\"\t\"I/O-Zugriffe\"\t\"wahlfreie I/Os\"\t\"I/O-Volumen\"\t\"Rekursionstiefe\"\t\"Offene Dateien\"\t\"Idle-Zeit vor I/O\"\t\"Idle-Zeit nach I/O\"\n");
        }

        inline void dumpTime(_proFloat seconds) {
			if (seconds < 0) {
				fputc('-', out);
				seconds = -seconds;
			}
            int secs    = (int)seconds;
            int mins    = secs/60;  secs -= 60*mins;
            int hours   = mins/60;  mins -= 60*hours;
            fprintf(out, "%d:%02d:%02d", hours, mins, secs);
        }

        inline void dumpTimeEx(_proFloat seconds) {
            int milli   = (int)(seconds * 1000.0);
            int secs    = (int)seconds;
            int mins    = secs/60;  secs -= 60*mins;
            int hours   = mins/60;  mins -= 60*hours;
            fprintf(out, "%d:%02d:%02d.%03d", hours, mins, secs, milli);
        }

        inline void dumpValue(_proTStates &stat, int valNum) {
			_proFloat f = stat[valNum];
            if ((_proValueType[valNum] & SEQAN_PROSTATE) == 0)
				f = _proData<>::_proValue[valNum] - f;

			switch (_proValueType[valNum] & SEQAN_PROTYPEMASK) {
				case SEQAN_PROINT:   									// state value -> print last seen maximum
					fprintf(out, "%.0f", f);
					break;

				case SEQAN_PROFLOAT:
					fprintf(out, "%f", f);
					break;

				case SEQAN_PROTIME:
					dumpTimeEx(f);
			}
        }

        inline void dumpSysValues(_proTStates &stat) {
            for(int i = 0; i < SEQAN_PROINDEXCOUNT - SEQAN_PROEXTRACOUNT; ++i) {
                dumpTab();
                dumpValue(stat, i);
            }
        }

        inline void dumpExtraValues(_proTStates &stat) {
            for(int i = 0; i < _proData<>::_proExtraCount; ++i) {
                dumpTab();
                dumpValue(stat, SEQAN_PROINDEXCOUNT - SEQAN_PROEXTRACOUNT + i);
            }
	}
	
        inline void dumpMark() {
            if (!mark.empty()) {
                dumpTab();
                fprintf(out, "\"%s\"", mark.c_str());
                mark.erase();
            }
        }

        inline void dump(_proTStates &stat) {
			setTime(_proData<>::_proValue);
            dumpNext += dumpStep;
            bol = true;
            bool _flush = ((dumpStep == 0.0)) || ((lines & 16) == 0);

            dumpSysValues(stat);
            dumpExtraValues(stat);
            dumpMark();
            dumpEndl();
            if (_flush) fflush(out);
            ++lines;
        }

        inline void signalDumpTest(_proFloat now) {
            if (dumpStep > 0 && now > dumpNext && running) {
                dump(last);
                maximize(all, last);
                sync(last);
            }
        }

        inline void signalNewMax(int valNum) {
            if (running)
                if (last[valNum] < _proData<>::_proValue[valNum])
                    last[valNum] = _proData<>::_proValue[valNum];
        }

        inline void setMark(const char *text) {
            if (running) {
                mark = text;
                if (dumpStep == 0.0) {
                    dump(last);                 // manual dump;
                    maximize(all, last);
                    sync(last);
                }
            }
        }
        
        inline void reset() {
        	syncTime(last);
        }

        inline void setEndMark(const char *text) {
            if (running) {
				setMark(text);
				reset();
			}
        }

    private:
        
        bool bol;   // begin of line
    };



/*
    inline void _proSignalDumpTest(_proFloat now);
    inline void _proSignalNewMax(int valNum);
    inline void _proMark(const char *text);
    inline void _proEndMark(const char *text);
    inline void _proReset();

    inline void _proSet(int valNum, _proFloat value);
    inline void _proAdd(int valNum, _proFloat value);
    inline void _proSub(int valNum, _proFloat value);
    
    // simple interface for external programs
    inline void *_proMalloc(size_t size);
    inline void _proFree(void *_ptr);
*/

    inline void _proSignalDumpTest(_proFloat now) {
        if (_proData<>::_proPFileStream) _proData<>::_proPFileStream->signalDumpTest(now);
    }

    inline void _proSignalNewMax(int valNum) {
        if (((_proValueType[valNum] & SEQAN_PROSTATE) != 0)) {
            if (_proData<>::_proPFileStream) _proData<>::_proPFileStream->signalNewMax(valNum);
            if (_proData<>::_proPFile)       _proData<>::_proPFile->signalNewMax(valNum);
        }
    }

    inline void _proMark(const char *text) {
        if (_proData<>::_proPFileStream) _proData<>::_proPFileStream->setMark(text);
        if (_proData<>::_proPFile)       _proData<>::_proPFile->setMark(text);
    }

    inline void _proEndMark(const char *text) {
        if (_proData<>::_proPFileStream) { _proData<>::_proPFileStream->setEndMark(text); }
        if (_proData<>::_proPFile)       { _proData<>::_proPFile->setEndMark(text); }
    }

    inline void _proReset() {
        if (_proData<>::_proPFileStream) { _proData<>::_proPFileStream->reset(); }
        if (_proData<>::_proPFile)       { _proData<>::_proPFile->reset(); }
    }




    inline void _proSet(_proValueIndex valNum, _proFloat value) {
        _proFloat now = sysTime();
        _proData<>::_proLastUpdate[valNum] = now;
        if (_proData<>::_proValue[valNum] < value) {
            _proData<>::_proValue[valNum] = value;
            _proSignalNewMax(valNum);
        } else
            _proData<>::_proValue[valNum] = value;
        _proSignalDumpTest(now);
    }

    inline void _proAdd(_proValueIndex valNum, _proFloat value) {
        _proFloat now = sysTime();
        _proData<>::_proValue[valNum] += value;
        _proData<>::_proLastUpdate[valNum] = now;
        if (valNum == SEQAN_PROIO) _proAdd(SEQAN_PROIORANDOM, 1);
        _proSignalNewMax(valNum);
        _proSignalDumpTest(now);
    }

    inline void _proSub(_proValueIndex valNum, _proFloat value) {
        _proFloat now = sysTime();
        _proData<>::_proValue[valNum] -= value;
        _proData<>::_proLastUpdate[valNum] = now;
        _proSignalDumpTest(now);
    }
    
    // simple interface for external programs
    inline void *_proMalloc(size_t size) {
    	size_t *ptr = reinterpret_cast<size_t*>(malloc(size + sizeof(size_t)));
    	if (ptr) {
    		_proAdd(SEQAN_PROMEMORY, *ptr = size);
//			printf("_proMalloc %x size %d\n", ptr, size);
    		++ptr;
    	}
    	return ptr;
    }

    inline void _proFree(void *_ptr) {
    	size_t *ptr = reinterpret_cast<size_t*>(_ptr);
    	if (ptr) {
    		--ptr;
//			printf("_proFree   %x size %d\n", _ptr, *ptr);
    		_proSub(SEQAN_PROMEMORY, *ptr);
    	}
    	free(ptr);
    }

#endif
