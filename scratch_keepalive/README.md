# scratch_keepalive

Keeps the Quest **b1042 genomics-scratch** tree alive against the 30-day
inactivity auto-purge (which deletes files not *accessed* in 30 days).

## How it works
`keepalive_touch.sh` is a self-resubmitting Slurm job:
1. Runs `find /gpfs/projects/b1042/epifluidlab/yoshii/ -exec touch -a -m {} +`
   to refresh access+modify times on the whole tree (~83k files).
2. Resubmits itself with `sbatch --begin=now+7days`, so it runs every 7 days
   forever — surviving login-node rotation (unlike `crontab`, and `scrontab`
   is disabled on Quest).

Account `b1042`, partition `genomics`, ~2 min runtime, 4G mem, 2h wall limit.

## Why a self-resubmitting Slurm job
- `scrontab` (Slurm cron) is **disabled on this cluster**.
- plain `crontab` is per-login-node; Quest rotates quser40-44, so it's unreliable.

## Operate it
- Check status / next run:  `squeue -u jmj7858 -n keepalive_touch`
- View last run log:        `ls -t logs/ | head; cat logs/<newest>.out`
- **Stop the chain:**       `scancel` the pending `keepalive_touch` job.
- Restart the chain:        `sbatch keepalive_touch.sh`

## ⚠️ This is a stopgap, not a backup
Files stay on unbacked-up scratch. For anything precious, move it to permanent
storage. See FILES.md / JOBS.md and the session note in ../.
