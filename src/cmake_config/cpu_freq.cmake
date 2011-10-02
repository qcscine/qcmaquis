# TODO: find a way to get this information during config
# Mac: sysctl -a hw
# Unix: cat /proc/cpuinfo
set(CPU_FREQ 2.8e9 CACHE STRING "CPU Frequency")
add_definitions(-DCPU_FREQ=${CPU_FREQ})
