#!/usr/bin/env bash
{
python scripts/plot.py --method fineeta --samples outer &
python scripts/plot.py --method fineeta --samples inner &
python scripts/plot.py --method ed --samples outer &
python scripts/plot.py --method ed --samples inner &

#Warn the user after all the jobs have finished
pids=()
for i in `jobs -p`; do
    pids+=("${i}")
done
for pid in ${pids[@]}; do
    wait ${pid}
done
echo "All jobs have finished."
}
exit $? #protects against run-time appends
