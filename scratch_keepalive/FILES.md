# FILES.md — scratch_keepalive/

| File | Purpose | How generated |
|------|---------|---------------|
| `keepalive_touch.sh` | Self-resubmitting Slurm job: touches `yoshii/` tree then re-submits itself +7 days. | Hand-written 2026-05-30. |
| `logs/keepalive_<jobid>.out` | Per-run stdout/stderr: timestamp, touch exit code, file count, resubmit result. | Emitted by each Slurm run. |
| `README.md` | What this folder does and how to operate/stop it. | Hand-written 2026-05-30. |
| `JOBS.md` | Submitted-job tracking (job id, command, log path, next run). | Hand-written 2026-05-30. |
| `FILES.md` | This file. | Hand-written 2026-05-30. |
