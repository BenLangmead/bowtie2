/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SPINLOCK_H_
#define SPINLOCK_H_

/**
 * This non-reentrant spinlock implementation for i386 and x86_64 is
 * based on free code by Gert Boddaert:
 *
 * http://www.codeproject.com/KB/threads/spinlocks.aspx
 *
 * Using spinlocks instead of the heavier pthreads mutexes can, in some
 * cases, help Bowtie perform better for large numbers of threads.
 */

// (If the user hasn't specified this, then there's no need for any
// kind of locking, so we skip this header)
#ifdef BOWTIE_PTHREADS

#if defined(__GNUC__)
#if defined(__x86_64__) || defined(__i386__)
#ifdef TRY_SPINLOCK
#define USE_SPINLOCK
#endif
#endif
#endif

#ifdef USE_SPINLOCK

#if defined(__x86_64__)
#define SPINLOCK_WORD long
#else
#define SPINLOCK_WORD int
#endif

class SpinLock {
  public:    // inlined constructor

    // inlined NON-virtual destructor
    inline SpinLock() : m_s(1) {}
    inline ~SpinLock() {}

    // enter the lock, spinlocks (with/without Sleep)
    // when mutex is already locked
    inline void Enter(void)
    {
    	SPINLOCK_WORD prev_s;
        do
        {
            prev_s = TestAndSet(&m_s, 0);
            if (m_s == 0 && prev_s == 1)
            {
            	// The lock and was unlocked and we grabbed it
                break;
            }
            // reluinquish current timeslice (can only
            // be used when OS available and
            // we do NOT want to 'spin')
            // HWSleep(0);
        }
        while (true);
    }
    // Tries to enter the lock, returns 0
    // when mutex is already locked,
    // returns != 0 when success
    inline int TryEnter(void)
    {
    	SPINLOCK_WORD prev_s = TestAndSet(&m_s, 0);
        if (m_s == 0 && prev_s == 1)
        {
            return 1;
        }
        return 0;
    }

    // Leaves or unlocks the mutex
    // (should only be called by lock owner)
    inline void Leave(void)
    {
        TestAndSet(&m_s, 1);
    }

  protected:
    // sets BIT value and returns previous
    // value.in 1 atomic un-interruptable operation
	  SPINLOCK_WORD TestAndSet(SPINLOCK_WORD* pTargetAddress, SPINLOCK_WORD nValue);

  private:
	SPINLOCK_WORD m_s;
};

// This part is Platform dependent!

/* The following piece of code can be found
 in function AtomicExchange of file atomicops-internals-x86.h
 under http://google-perftools.googlecode.com/svn/trunk/src/base/
*/
#if defined(__x86_64__)
#define TAS(_lw, _res) \
 asm volatile("xchgq %1,%0":"=r"(_res):"m"(*_lw),"0"(_res):"memory")
#elif defined(__i386__)
#define TAS(_lw, _res) \
 asm volatile("xchgl %1,%0":"=r"(_res):"m"(*_lw),"0"(_res):"memory")
#else
#error "Architecture is neither x86_64 nor i386 (see spinlock.h)"
#endif
/* TAS is only defined for GNUC on x86_64 and i386,
 that is where GNUC x86 inline assembly can be used */

inline SPINLOCK_WORD SpinLock::TestAndSet(SPINLOCK_WORD* pTargetAddress, SPINLOCK_WORD nValue) {

#if 0
    __asm
    {
        mov edx, dword ptr [pTargetAddress]
        mov eax, nValue
        lock xchg eax, dword ptr [edx]
    }
#else
    TAS(pTargetAddress, nValue);
    return nValue;
#endif
    // mov = 1 CPU cycle
    // lock = 1 CPU cycle
    // xchg = 3 CPU cycles
}
#endif /*USE_SPINLOCK*/

#endif /*BOWTIE_PTHREADS*/

#endif /*SPINLOCK_H_*/
