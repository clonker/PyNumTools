#!/usr/bin/env bash

set +e

nosetests -s -vv --with-doctest --doctest-options=+NORMALIZE_WHITESPACE,+ELLIPSIS
