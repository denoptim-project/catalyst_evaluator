# Useful commands

Here we collect commands that are useful in one way or another when performing the catalyst evaluations.

* Count the jobs that are waiting for the termination of something (run this from the folder from which all such jobs have been sumbitted to the catalysts evaluator):
```
echo jobs_waiting=$(tail -n 1  *-*/*FProv*.log | grep -c 'Still waiting for')
```

* List the active processes that own a locked seat for conformational search 
```
ls /tmp/activeConfSearchOrRefine* | while read f ; do p=$(cat $f | sed 's/:/ /' | sed 's/)/ /' | awk '{print $3}'); echo $f $p ;  ps f -p $p ; done
```



