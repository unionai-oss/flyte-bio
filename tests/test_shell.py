"""Tests for shell-task behaviour the modules rely on."""

import flyte
from flyte.extras import shell
from flyte.io import File

from flyte_bio.limits import JAVA_HEAP_MB, TASK_LIMITS
from flyte_bio.modules.gunzip import GUNZIP_IMAGE
from tests.framework import env

# Writes the container's CPU and memory limits as cgroup v2 reports them
# ("<quota> <period>" and bytes), or cgroup v1's equivalents.
limits_cmd = shell.create(
    name="resource_limits",
    image=GUNZIP_IMAGE,
    resources=flyte.Resources(cpu=1, memory="1Gi"),
    outputs={"limits": File},
    script=r"""
        if [ -f /sys/fs/cgroup/cpu.max ]; then
            CPU=$(cat /sys/fs/cgroup/cpu.max)
            MEM=$(cat /sys/fs/cgroup/memory.max)
        else
            CPU="$(cat /sys/fs/cgroup/cpu/cpu.cfs_quota_us) $(cat /sys/fs/cgroup/cpu/cpu.cfs_period_us)"
            MEM=$(cat /sys/fs/cgroup/memory/memory.limit_in_bytes)
        fi
        printf '%s\n%s\n' "$CPU" "$MEM" > {outputs.limits}
    """,
)

# Writes TASK_LIMITS' CPUs and memory, and the JVM heap they give.
task_limits_cmd = shell.create(
    name="task_limits",
    image=GUNZIP_IMAGE,
    resources=flyte.Resources(cpu=1, memory="1Gi"),
    # Resources aren't part of the cache key: the overridden call would reuse the first one's output.
    cache="disable",
    outputs={"limits": File},
    script=rf"""
        {TASK_LIMITS}
        echo "$TASK_CPUS $TASK_MEMORY_MB {JAVA_HEAP_MB}" > {{outputs.limits}}
    """,
)

env.add_dependency(flyte.TaskEnvironment.from_task("resource_limits", limits_cmd.as_task()))
env.add_dependency(flyte.TaskEnvironment.from_task("task_limits", task_limits_cmd.as_task()))


async def read_text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


@env.task
async def test_override_resources() -> None:
    task = limits_cmd.override(resources=flyte.Resources(cpu=2, memory="3Gi"))
    cpu, mem = (await read_text(await task())).split("\n")[:2]
    quota, period = (int(x) for x in cpu.split())
    assert quota / period == 2, f"CPU limit {cpu!r}, expected 2 CPUs"
    assert int(mem) == 3 * 1024**3, f"memory limit {mem!r}, expected 3Gi"


@env.task
async def test_task_limits() -> None:
    # A tool's threads and memory follow the resources the task got, overridden or not.
    limits = (await read_text(await task_limits_cmd())).split()
    assert limits == ["1", "1024", "819"], limits
    task = task_limits_cmd.override(resources=flyte.Resources(cpu=1.5, memory="3Gi"))
    limits = (await read_text(await task())).split()
    assert limits == ["2", "3072", "2457"], limits


tests = [test_override_resources, test_task_limits]
