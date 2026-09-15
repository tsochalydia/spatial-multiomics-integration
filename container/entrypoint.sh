#!/usr/bin/env bash
set -euo pipefail

if [ -z "${HOME:-}" ] || [ "$HOME" = "/" ] || [ ! -w "${HOME:-/nonexistent}" ]; then
    export HOME="/tmp/home-$(id -u)"
    mkdir -p "$HOME"
fi

eval "$(micromamba shell hook --shell bash)"

exec "$@"
