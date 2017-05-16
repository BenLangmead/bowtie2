import sys

# gets the average time from this file where
# each line is a different time
# taken from usr/bin/time -v.
def get_avg_time(timeFile, run):
    totalA_min = 0
    totalA_hr = 0
    totalA_sec = 0.0
    for line in timeFile:
        data = line.strip().split()
        time = data[-1]
        if time.count(':') > 1:
            hrs = int(time[0])
            mins = int(time[time.index(':') + 1 : time.index(':') + 3])
            secs = float(time[time.index(':') + 4:])
        else:
            hrs = 0
            mins = int(time[0:time.index(':')])
            secs = float(time[time.index(':') + 1 : ])
            totalA_min += mins
            totalA_hr += hrs
            totalA_sec += secs
            
    totalA_min += totalA_hr * 60
    totalA_sec += totalA_min * 60
    
    avgA_time = totalA_sec / 3.0
    Atime = int(avgA_time)
    A_sec_dec = avgA_time - Atime
    A_sec = Atime % 60
    Atime //= 60
    A_min = Atime % 60
    A_hr = Atime // 60
    
    print(run)
    print('\nElapsed (wall clock) time (h:mm:ss): {:02d}:{:02d}:{:05.2f}\n'.format(A_hr, A_min, A_sec + A_sec_dec))    

# gets the average memory footprint from this file where
# each line is a different memory footprint
# taken from usr/bin/time -v.
def get_avg_mem(memFile, run):
    total_mem = 0
    for line in memFile:
        data = line.strip().split()
        total_mem += int(data[-1])

    total_mem //= 3
    print(run)
    print('\nMaximum resident set size (kbytes): {}\n'.format(total_mem))    

# takes 4 command line inputs: 1 = time file for run A,
# 2 = memory file for run A, 3 = time file for run B,
# 4 = memory file for run B.
def main():
    runA_times = open(sys.argv[1], "r")
    runA_mem = open(sys.argv[2], "r")
    runB_times = open(sys.argv[3], "r")
    runB_mem = open(sys.argv[4], "r")
    
    print('\n\n\n\n')
    get_avg_time(runA_times, 'run A time:')
    get_avg_time(runB_times, 'run B time:')
    get_avg_mem(runA_mem, 'run A peak memory footprint:')
    get_avg_mem(runB_mem, 'run B peak memory footprint:')
    
main()
