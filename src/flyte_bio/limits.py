"""The container's CPU and memory limits, read when a shell task runs.

A tool's thread count and memory (``-t``, ``-Xmx``, sort buffers) must follow
the resources the task actually got, not the module's defaults, so that a
``task.override(resources=...)`` reaches the tool too. Upstream's modules get
the same from Nextflow's ``task.cpus`` and ``task.memory``.

Start a shell task's script with :data:`TASK_LIMITS` to set two variables:

- ``TASK_CPUS``: the CPU limit, rounded up (at least 1).
- ``TASK_MEMORY_MB``: the memory limit in MiB.

They come from the cgroup (v2, else v1). Without a limit they fall back to the
machine's CPUs and memory.
"""

TASK_LIMITS = r"""
        if [ -r /sys/fs/cgroup/cpu.max ]; then
            read -r CPU_QUOTA CPU_PERIOD < /sys/fs/cgroup/cpu.max
            MEMORY_LIMIT=$(cat /sys/fs/cgroup/memory.max)
        else
            CPU_QUOTA=$(cat /sys/fs/cgroup/cpu/cpu.cfs_quota_us 2>/dev/null || echo max)
            CPU_PERIOD=$(cat /sys/fs/cgroup/cpu/cpu.cfs_period_us 2>/dev/null || echo 100000)
            MEMORY_LIMIT=$(cat /sys/fs/cgroup/memory/memory.limit_in_bytes 2>/dev/null || echo max)
        fi
        if [ "$CPU_QUOTA" = max ] || [ "$CPU_QUOTA" -le 0 ]; then
            TASK_CPUS=$(nproc 2>/dev/null || echo 1)
        else
            TASK_CPUS=$(( (CPU_QUOTA + CPU_PERIOD - 1) / CPU_PERIOD ))
        fi
        HOST_MEMORY_MB=0
        while read -r KEY VALUE UNIT; do
            if [ "$KEY" = MemTotal: ]; then HOST_MEMORY_MB=$(( VALUE / 1024 )); fi
        done < /proc/meminfo
        TASK_MEMORY_MB=$HOST_MEMORY_MB
        if [ "$MEMORY_LIMIT" != max ]; then
            LIMIT_MB=$(( MEMORY_LIMIT / 1048576 ))
            if [ "$HOST_MEMORY_MB" -eq 0 ] || [ "$LIMIT_MB" -lt "$HOST_MEMORY_MB" ]; then TASK_MEMORY_MB=$LIMIT_MB; fi
        fi
""".strip()

# As upstream: a JVM tool's heap (-Xmx) is 80% of the task's memory, in MiB.
# Use it in a script that starts with TASK_LIMITS.
JAVA_HEAP_MB = "$(( TASK_MEMORY_MB * 8 / 10 ))"
