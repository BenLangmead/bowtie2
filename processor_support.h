#ifndef PROCESSOR_SUPPORT_H_
#define PROCESSOR_SUPPORT_H_

// Utility class ProcessorSupport provides POPCNTenabled() to determine
// processor support for POPCNT instruction. It uses CPUID to
// retrieve the processor capabilities.
// for Intel ICC compiler __cpuid() is an intrinsic 
// for Microsoft compiler __cpuid() is provided by #include <intrin.h>
// for GCC compiler __get_cpuid() is provided by #include <cpuid.h>

// Intel compiler defines __GNUC__, so this is needed to disambiguate

#if defined(__INTEL_COMPILER)
#   define USING_INTEL_COMPILER
#elif defined(__GNUC__)
#   define USING_GCC_COMPILER
#   include <cpuid.h>
#elif defined(_MSC_VER)
// __MSC_VER defined by Microsoft compiler
#define USING MSC_COMPILER
#endif

struct regs_t {unsigned int EAX, EBX, ECX, EDX;};
#define BIT(n) ((1<<n))

class ProcessorSupport {

#ifdef POPCNT_CAPABILITY 

public: 
    ProcessorSupport() { } 
    bool POPCNTenabled()
    {
    // from: Intel® 64 and IA-32 Architectures Software Developer’s Manual, 325462-036US,March 2013
    //Before an application attempts to use the POPCNT instruction, it must check that the
    //processor supports SSE4.2
    //“(if CPUID.01H:ECX.SSE4_2[bit 20] = 1) and POPCNT (if CPUID.01H:ECX.POPCNT[bit 23] = 1)”
    //
    // see p.272 of http://download.intel.com/products/processor/manual/253667.pdf available at
    // http://www.intel.com/content/www/us/en/processors/architectures-software-developer-manuals.html
    // Also http://en.wikipedia.org/wiki/SSE4 talks about available on Intel & AMD processors

    regs_t regs;

    try {
#if ( defined(USING_INTEL_COMPILER) || defined(USING_MSC_COMPILER) )
        __cpuid((int *) &regs,0); // test if __cpuid() works, if not catch the exception
        __cpuid((int *) &regs,0x1); // POPCNT bit is bit 23 in ECX
#elif defined(USING_GCC_COMPILER)
        __get_cpuid(0x1, &regs.EAX, &regs.EBX, &regs.ECX, &regs.EDX);
#else
        std::cerr << “ERROR: please define __cpuid() for this build.\n”; 
        assert(0);
#endif
        if( !( (regs.ECX & BIT(20)) && (regs.ECX & BIT(23)) ) ) return false;
    }
    catch (int e) {
        return false;
    }
    return true;
    }

#endif // POPCNT_CAPABILITY
};

#endif /*PROCESSOR_SUPPORT_H_*/




