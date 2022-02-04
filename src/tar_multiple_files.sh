#!/bin/bash

for f in "$@"
do
  sudo tar "$f".tar "$f"
done