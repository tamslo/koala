# Commands for manual execution

Docker fails to execute some commands with the error `Error grabbing logs: unexpected EOF`.
This happens for NovoAlign and BEDTools `genomecov`.

So we are running the commands manually in the regarding docker container and write stats to a log file (`docker stats <container id> > <stats_file>`)
that are evaluated with the `docker_stats_log_to_runstats.py` script.

With SSH connections, it is useful to run commands as jobs (append a `&`) and write a log file (`command > logfile 2>&1`) so that the command does not fail when the connection is lost.
