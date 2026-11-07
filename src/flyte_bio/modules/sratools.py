"""sratools — download runs from NCBI's Sequence Read Archive and dump them to FASTQ with SRA Toolkit.

Exposes three functions:

- :func:`sratools_ncbi_settings` writes the SRA Toolkit's user settings file
  (``*.mkfg``), which :func:`sratools_prefetch` and
  :func:`sratools_fasterqdump` both need (``NCBI_SETTINGS``). Without
  ``existing_settings`` it generates ``user-settings.mkfg`` with a fresh
  ``/LIBS/GUID`` and ``/libs/cloud/report_instance_identity = "true"``, as
  upstream does when no settings exist. With ``existing_settings`` (a
  settings file of your own, what upstream finds through ``NCBI_SETTINGS``
  or ``vdb-config``), it checks that
  ``prefetch`` and ``fasterq-dump`` accept it and that it is named
  ``*.mkfg``, and returns it under its own name. The task is cached, so the
  generated GUID stays the same across runs.
- :func:`sratools_prefetch` downloads one accession (``SRR...``, ``DRR...``,
  ``ERR...``) into a Dir named after it, holding ``<id>.sra`` or
  ``<id>.sralite``. As upstream, prefetch is retried with exponential
  backoff (``retry_args``: ``<attempts> <base delay s> <max delay s>``,
  default ``5 1 100``), and the download is checked with ``vdb-validate``,
  or, when the archive has no checksums, against NCBI's md5s.
- :func:`sratools_fasterqdump` dumps a :func:`sratools_prefetch` Dir to
  gzipped FASTQ: ``<prefix>.fastq.gz`` (``single_end``) or
  ``<prefix>_1.fastq.gz`` / ``<prefix>_2.fastq.gz``. For paired-end data the
  unpaired reads fasterq-dump writes to ``<prefix>.fastq`` are dropped, as
  upstream. ``args`` goes to fasterq-dump and ``args2`` to pigz.

``certificate`` (prefetch and fasterq-dump) gives access to protected data: a
``.jwt`` cart file is passed as ``--perm``, a ``.ngc`` repository key as
``--ngc``. fasterq-dump and pigz use the task's CPUs.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

SRATOOLS_IMAGE = "community.wave.seqera.io/library/sra-tools:3.2.1--c221ea63b92c74bf"
# prefetch's image adds curl, to fetch NCBI's md5s; fasterq-dump's adds pigz.
SRATOOLS_CURL_IMAGE = "community.wave.seqera.io/library/sra-tools_curl:c1e36977de327cd2"
SRATOOLS_PIGZ_IMAGE = "community.wave.seqera.io/library/sra-tools_pigz:dd7a0b6dd85ec25a"
SETTINGS_RESOURCES = flyte.Resources(cpu=1, memory="1Gi")
# Upstream labels prefetch process_low and fasterq-dump process_medium. Large
# runs need more disk: fasterq-dump's temporary files are several times the
# size of the .sra.
PREFETCH_RESOURCES = flyte.Resources(cpu=1, memory="2Gi", disk="20Gi")
FASTERQDUMP_RESOURCES = flyte.Resources(cpu=4, memory="8Gi", disk="50Gi")
PREFETCH_RETRY_ARGS = "5 1 100"

# The certificate flag for a .jwt cart file or .ngc repository key in CERT[0], as upstream.
CERTIFICATE_ARGS = r"""
        CERT_ARGS=()
        if [ ${#CERT[@]} -gt 0 ]; then
            case "${CERT[0]}" in
                *.jwt) CERT_ARGS=(--perm "${CERT[0]}") ;;
                *.ngc) CERT_ARGS=(--ngc "${CERT[0]}") ;;
            esac
        fi
""".strip()

sratools_ncbi_settings_cmd = shell.create(
    name="sratools_ncbi_settings",
    image=SRATOOLS_IMAGE,
    resources=SETTINGS_RESOURCES,
    inputs={"existing_settings": File | None},
    outputs={"results": Dir},
    script=r"""
        shopt -s nullglob; EXISTING=({inputs.existing_settings}); shopt -u nullglob
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"

        UUID=$(cat /proc/sys/kernel/random/uuid)
        CONFIG="/LIBS/GUID = \"$UUID\"\n/libs/cloud/report_instance_identity = \"true\"\n"

        # Without settings, create a file suitable for prefetch and fasterq-dump.
        # Existing settings that lack the required values are an error. Upstream
        # finds existing settings with vdb-config, but in a fresh container that
        # call itself writes a settings file with only a GUID, which would then
        # be taken for existing settings; so only existing_settings counts here.
        if [ ${#EXISTING[@]} -eq 0 ]; then
            printf "$CONFIG" > 'user-settings.mkfg'
        else
            export NCBI_SETTINGS="${EXISTING[0]}"
            prefetch --help &> /dev/null
            if [[ $? = 78 ]]; then
                echo "You have an existing vdb-config at '${NCBI_SETTINGS}' but it is"\
                "missing the required entries for /LIBS/GUID and"\
                "/libs/cloud/report_instance_identity."\
                "Feel free to add the following to your settings file:" >&2
                echo "$(printf "$CONFIG")" >&2
                exit 1
            fi
            fasterq-dump --help &> /dev/null
            if [[ $? = 78 ]]; then
                echo "You have an existing vdb-config at '${NCBI_SETTINGS}' but it is"\
                "missing the required entries for /LIBS/GUID and"\
                "/libs/cloud/report_instance_identity."\
                "Feel free to add the following to your settings file:" >&2
                echo "$(printf "$CONFIG")" >&2
                exit 1
            fi
            if [[ "${NCBI_SETTINGS}" != *.mkfg ]]; then
                echo "The detected settings '${NCBI_SETTINGS}' do not have the required"\
                "file extension '.mkfg'." >&2
                exit 1
            fi
            cp "${NCBI_SETTINGS}" ./
        fi
        mv ./*.mkfg "$OUT"/
    """,
)

sratools_prefetch_cmd = shell.create(
    name="sratools_prefetch",
    image=SRATOOLS_CURL_IMAGE,
    resources=PREFETCH_RESOURCES,
    inputs={"id": str, "ncbi_settings": File, "certificate": File | None, "retry_args": str, "args": str},
    defaults={"retry_args": PREFETCH_RETRY_ARGS, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        SETTINGS=({{inputs.ncbi_settings}})
        shopt -s nullglob; CERT=({{inputs.certificate}}); shopt -u nullglob
        ID={{inputs.id}}
        RETRY_ARGS={{inputs.retry_args}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        {CERTIFICATE_ARGS}

        retry_with_backoff() {{
            local max_attempts=${{1}}
            local delay=${{2}}
            local max_time=${{3}}
            local attempt=1
            local output=
            local status=

            # Remove the first three arguments to this function in order to access
            # the 'real' command with `${{@}}`.
            shift 3

            while [ ${{attempt}} -le ${{max_attempts}} ]; do
                output=$("${{@}}")
                status=${{?}}

                if [ ${{status}} -eq 0 ]; then
                    break
                fi

                if [ ${{attempt}} -lt ${{max_attempts}} ]; then
                    echo "Failed attempt ${{attempt}} of ${{max_attempts}}. Retrying in ${{delay}} s." >&2
                    sleep ${{delay}}
                elif [ ${{attempt}} -eq ${{max_attempts}} ]; then
                    echo "Failed after ${{attempt}} attempts." >&2
                    return ${{status}}
                fi

                attempt=$(( ${{attempt}} + 1 ))
                delay=$(( ${{delay}} * 2 ))
                if [ ${{delay}} -ge ${{max_time}} ]; then
                    delay=${{max_time}}
                fi
            done

            echo "${{output}}"
        }}

        cp "${{SETTINGS[0]}}" .
        export NCBI_SETTINGS="$PWD/$(basename "${{SETTINGS[0]}}")"

        retry_with_backoff $RETRY_ARGS \
            prefetch \
            $ARGS "${{CERT_ARGS[@]}}" \
            "$ID" \
            || exit 1

        # Check file integrity using vdb-validate or (when the archive contains no checksums) md5sum.
        vdb-validate "$ID" > vdb-validate_result.txt 2>&1 || exit 1
        if grep -q "checksums missing" vdb-validate_result.txt; then
            VALID_MD5SUMS=$(curl --silent --fail --location --retry 3 --retry-delay 60 \
                "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?filetype=run&acc=$ID")
            LOCAL_MD5SUMS=$(md5sum "$ID"/* | cut -f1 -d' ')
            if ! grep -q -F -f <(echo "$LOCAL_MD5SUMS") <(echo "$VALID_MD5SUMS"); then
                echo "MD5 sum check failed" 1>&2
                exit 1
            fi
        fi
        mv "$ID" "$OUT"/ || exit 1
    """,
)

sratools_fasterqdump_cmd = shell.create(
    name="sratools_fasterqdump",
    image=SRATOOLS_PIGZ_IMAGE,
    resources=FASTERQDUMP_RESOURCES,
    inputs={
        "sra": Dir,
        "ncbi_settings": File,
        "certificate": File | None,
        "prefix": str,
        "single_end": bool,
        "args": str,
        "args2": str,
    },
    defaults={"args": "", "args2": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        SETTINGS=({{inputs.ncbi_settings}})
        shopt -s nullglob; CERT=({{inputs.certificate}}); shopt -u nullglob
        SRA={{inputs.sra}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        {CERTIFICATE_ARGS}

        cp "${{SETTINGS[0]}}" .
        export NCBI_SETTINGS="$PWD/$(basename "${{SETTINGS[0]}}")"

        # As upstream, fasterq-dump reads the accession's directory by its plain
        # name: a '.' in the path trips it (https://github.com/ncbi/sra-tools/issues/865).
        shopt -s nullglob; RUNS=("$SRA"/*.sra "$SRA"/*.sralite); shopt -u nullglob
        if [ ${{#RUNS[@]}} -eq 0 ]; then
            echo "no .sra or .sralite file in the input directory" >&2
            exit 1
        fi
        ACC=$(basename "${{RUNS[0]}}"); ACC=${{ACC%.sralite}}; ACC=${{ACC%.sra}}
        ln -s "$SRA" "$ACC"

        fasterq-dump \
            $ARGS \
            --threads $TASK_CPUS \
            --outfile "$P.fastq" \
            "${{CERT_ARGS[@]}}" \
            "$ACC" \
            || exit 1

        # Paired-end: drop the unpaired reads' file from the outputs.
        if [ {{inputs.single_end}} = false ]; then mv "$P.fastq" "$P" || echo 'No third file'; fi

        pigz \
            $ARGS2 \
            --no-name \
            --processes $TASK_CPUS \
            ./*.fastq \
            || exit 1
        mv ./*.fastq.gz "$OUT"/
    """,
)

sratools_env = flyte.TaskEnvironment.from_task("sratools_tools", sratools_ncbi_settings_cmd.as_task())
sratools_curl_env = flyte.TaskEnvironment.from_task("sratools_curl", sratools_prefetch_cmd.as_task())
sratools_pigz_env = flyte.TaskEnvironment.from_task("sratools_pigz", sratools_fasterqdump_cmd.as_task())
env = flyte.TaskEnvironment(
    name="sratools",
    depends_on=[sratools_env, sratools_curl_env, sratools_pigz_env],
)


async def sratools_ncbi_settings(existing_settings: File | None = None) -> File:
    """The SRA Toolkit settings file: generated, or ``existing_settings`` once checked."""
    results = await sratools_ncbi_settings_cmd(existing_settings=existing_settings)
    settings = [f async for f in results.walk() if f.name.endswith(".mkfg")]
    if len(settings) != 1:
        raise FileNotFoundError(f"sratools_ncbi_settings wrote {len(settings)} .mkfg files, expected 1")
    return settings[0]


async def sratools_prefetch(
    id: str,
    ncbi_settings: File,
    certificate: File | None = None,
    retry_args: str = PREFETCH_RETRY_ARGS,
    args: str = "",
) -> Dir:
    """Download the SRA accession ``id``: a Dir named ``id`` holding ``<id>.sra`` or ``<id>.sralite``."""
    results = await sratools_prefetch_cmd(
        id=id, ncbi_settings=ncbi_settings, certificate=certificate, retry_args=retry_args, args=args
    )
    return Dir.from_existing_remote(f"{results.path.rstrip('/')}/{id}")


async def sratools_fasterqdump(
    sra: Dir,
    ncbi_settings: File,
    prefix: str,
    single_end: bool,
    certificate: File | None = None,
    args: str = "",
    args2: str = "",
) -> list[File]:
    """Dump ``sra`` to ``<prefix>.fastq.gz`` (single-end) or ``<prefix>_1.fastq.gz`` and ``<prefix>_2.fastq.gz``."""
    results = await sratools_fasterqdump_cmd(
        sra=sra,
        ncbi_settings=ncbi_settings,
        certificate=certificate,
        prefix=prefix,
        single_end=single_end,
        args=args,
        args2=args2,
    )
    reads = [f async for f in results.walk() if f.name.endswith(".fastq.gz")]
    if not reads:
        raise FileNotFoundError("sratools_fasterqdump wrote no .fastq.gz")
    return sorted(reads, key=lambda f: f.name)
