# JOBS.md — scratch_keepalive/

Self-resubmitting weekly keep-alive chain (each run submits the next).

## Latest / next scheduled job
- **Next run job ID:** 9681601 — PENDING, begins **2026-06-06T21:47** (reason BeginTime)
- **Name:** keepalive_touch
- **Account/Partition:** b1042 / genomics
- **Log (when it runs):** `logs/keepalive_9681601.out`
- **Behavior:** touches `/gpfs/projects/b1042/epifluidlab/yoshii/`, then auto-submits the next run for +7 days (perpetual chain).

## Check status
```
squeue -u jmj7858 -n keepalive_touch
cat logs/keepalive_<jobid>.out
```

## Stop
`scancel` the pending `keepalive_touch` job (this breaks the chain; no further runs).

## History
| Date | Job ID | Result |
|------|--------|--------|
| 2026-05-30 | 9681588 | DONE: touched 83,637 files (exit 0), resubmitted 9681601 for 2026-06-06 |
| 2026-06-06 | 9681601 | PENDING (BeginTime) |
