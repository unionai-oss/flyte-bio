"""Tests for shell-task behaviour the modules rely on."""

import flyte
from flyte.extras import shell
from flyte.io import File

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

env.add_dependency(flyte.TaskEnvironment.from_task("resource_limits", limits_cmd.as_task()))


@env.task
async def test_override_resources() -> None:
    task = limits_cmd.as_task().override(resources=flyte.Resources(cpu=2, memory="3Gi"))
    limits = await task()
    async with limits.open("rb") as fh:
        cpu, mem = bytes(await fh.read()).decode().split("\n")[:2]
    quota, period = (int(x) for x in cpu.split())
    assert quota / period == 2, f"CPU limit {cpu!r}, expected 2 CPUs"
    assert int(mem) == 3 * 1024**3, f"memory limit {mem!r}, expected 3Gi"


tests = [test_override_resources]
