#!/usr/bin/env bash
set -e
# Usage:
#   check_version.sh
#
# Reads version from Cargo.toml and checks it against tags
#
# Credit to @richfitz for this from the dust package:
# https://github.com/mrc-ide/dust/blob/master/scripts/version_check
VERSION=${1:-$(grep '^version' Cargo.toml  | sed 's/.*= *//' | sed 's/"//g')}
TAG="v${VERSION}"

echo "Proposed version number '$VERSION'"

if echo "$VERSION" | grep -Eq "[0-9]+[.][0-9]+[.][0-9]+"; then
    echo "[OK] Version number in correct format"
else
    echo "[ERROR] Invalid format version number '$VERSION' must be in format 'x.y.z'"
    exit 1
fi

EXIT_CODE=0

echo "Updating remote git data"
git fetch --quiet

BRANCH_DEFAULT=$(git remote show origin | awk '/HEAD branch/ {print $NF}')
LAST_TAG=$(git describe --tags --abbrev=0 "origin/${BRANCH_DEFAULT}")

echo "Last tag was $LAST_TAG"

if git rev-parse "$TAG" >/dev/null 2>&1; then
    echo "[ERROR] Tag $TAG already exists - update version number in Cargo.toml"
    exit 1
else
    echo "[OK] Version number not yet present as git tag"
fi

MAJOR=$(echo $VERSION | cut -d. -f1)
MINOR=$(echo $VERSION | cut -d. -f2)
PATCH=$(echo $VERSION | cut -d. -f3)

LAST_VERSION=$(echo "$LAST_TAG" | sed 's/^v//')
LAST_MAJOR=$(echo $LAST_VERSION | cut -d. -f1)
LAST_MINOR=$(echo $LAST_VERSION | cut -d. -f2)
LAST_PATCH=$(echo $LAST_VERSION | cut -d. -f3)

if (( $MAJOR > $LAST_MAJOR )); then
    echo "[OK] Increasing MAJOR version"
    exit $EXIT_CODE
elif (( $MINOR > $LAST_MINOR )); then
    echo "[OK] Increasing MINOR version"
    exit $EXIT_CODE
elif (( $PATCH > $LAST_PATCH )); then
    echo "[OK] Increasing PATCH version"
    exit $EXIT_CODE
else
    echo "[ERROR] Version number has not increased relative to $LAST_VERSION"
    exit 1
fi
