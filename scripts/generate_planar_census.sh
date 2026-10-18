#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    printf 'Usage: %s OUTPUT_DIRECTORY\n' "$0" >&2
    exit 2
fi

for executable in plantri shortg mktemp sha256sum; do
    if ! command -v "$executable" >/dev/null 2>&1; then
        printf 'Required command not found: %s\n' "$executable" >&2
        exit 1
    fi
done

output_directory=$1
mkdir -p -- "$output_directory"
toolchain_file="$output_directory/toolchain.txt"
if [[ -e "$toolchain_file" ]]; then
    printf 'Refusing to overwrite: %s\n' "$toolchain_file" >&2
    exit 1
fi
temporary_file=
cleanup() {
    if [[ -n "$temporary_file" && -f "$temporary_file" ]]; then
        rm -- "$temporary_file"
    fi
}
trap cleanup EXIT

for n in 3 4 5 6; do
    destination="$output_directory/planar${n}.g6"
    if [[ -e "$destination" ]]; then
        printf 'Refusing to overwrite: %s\n' "$destination" >&2
        exit 1
    fi

    temporary_file=$(mktemp "$output_directory/.planar${n}.XXXXXX")
    plantri -p -c1 -g "$n" | shortg > "$temporary_file"
    line_count=$(wc -l < "$temporary_file")
    case "$n" in
        3) expected_count=2 ;;
        4) expected_count=6 ;;
        5) expected_count=20 ;;
        6) expected_count=99 ;;
    esac
    if [[ "$line_count" -ne "$expected_count" ]]; then
        printf 'Unexpected n=%s planar count: got %s, expected %s\n' \
            "$n" "$line_count" "$expected_count" >&2
        exit 1
    fi

    mv -- "$temporary_file" "$destination"
    temporary_file=
    printf 'n=%s: %s planar graph6 classes -> %s\n' \
        "$n" "$line_count" "$destination"
done

{
    printf 'plantri_command=plantri -p -c1 -g N\n'
    printf 'nauty_command=shortg\n'
    printf 'plantri_binary='; sha256sum "$(command -v plantri)"
    printf 'shortg_binary='; sha256sum "$(command -v shortg)"
    printf 'planar3='; sha256sum "$output_directory/planar3.g6"
    printf 'planar4='; sha256sum "$output_directory/planar4.g6"
    printf 'planar5='; sha256sum "$output_directory/planar5.g6"
    printf 'planar6='; sha256sum "$output_directory/planar6.g6"
} > "$toolchain_file"
printf 'Toolchain and output digests -> %s\n' "$toolchain_file"
