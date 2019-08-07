#!/usr/bin/env bash
python plotting/uproot_weights.py 3 inner &
python plotting/uproot_weights.py 4 inner &
python plotting/uproot_weights.py 5 inner &
python plotting/uproot_weights.py 6 inner &
python plotting/uproot_weights.py 3 outer &
python plotting/uproot_weights.py 4 outer &
python plotting/uproot_weights.py 5 outer &
python plotting/uproot_weights.py 6 outer &

python plotting/uproot_complete_fraction.py 3 inner &
python plotting/uproot_complete_fraction.py 4 inner &
python plotting/uproot_complete_fraction.py 5 inner &
python plotting/uproot_complete_fraction.py 6 inner &
python plotting/uproot_complete_fraction.py 3 outer &
python plotting/uproot_complete_fraction.py 4 outer &
python plotting/uproot_complete_fraction.py 5 outer &
python plotting/uproot_complete_fraction.py 6 outer &

python plotting/uproot_resolution.py 3 inner &
python plotting/uproot_resolution.py 4 inner &
python plotting/uproot_resolution.py 5 inner &
python plotting/uproot_resolution.py 6 inner &
python plotting/uproot_resolution.py 3 outer &
python plotting/uproot_resolution.py 4 outer &
python plotting/uproot_resolution.py 5 outer &
python plotting/uproot_resolution.py 6 outer &
