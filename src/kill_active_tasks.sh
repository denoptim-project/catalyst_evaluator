#!/bin/bash
#
# Kills all the tasks that hold lock files.
#
for actFile in $(ls /tmp/activeConfSearchOrRefine* /tmp/activeDFTTask* /tmp/activeXTBTask*  2> /dev/null)
do 
  pid="$(awk -F':|)' '{print $2}' "$actFile" )" 
  for spid in $(pgrep -P "$pid")
  do 
    kill -9 "$spid"
  done 
  lockFile="$(awk  '{print $5}' "$actFile")"
  rm -f "$lockFile"
  rm -f "$actFile"
  kill -9 "$pid" 
done
