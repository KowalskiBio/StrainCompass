# bactiment deployment guide (VM: proxmox1)

This VM is PRODUCTION and hosts other people's apps. Every deployment
follows this document. If any step would conflict with another app,
STOP and reconsider instead of forcing it.

## 0. Golden safety rules (non negotiable)

1. NEVER stop, restart, kill or "fix" any process that does not belong
   to bactiment. No `kill` on foreign PIDs, no `systemctl restart` of
   foreign units, no `pm2 restart`, nothing.
2. NEVER edit, overwrite or delete another app's files, nginx site
   configs, or data directories. bactiment gets its own files only.
3. NEVER bind a port that is already in use, and NEVER expose the app
   directly on 0.0.0.0. bactiment binds 127.0.0.1 only; the public side
   is the existing nginx.
4. ALWAYS back up bactiment's data before every deploy/push of a new
   version (see section 5). Users create real projects in it; their
   data outranks the release.
5. ALWAYS run the pre-deployment check (section 3) on the day of the
   deploy, even if you deployed yesterday. Things change.
6. When a command needs sudo, run it manually and consciously. The
   deploy user has no passwordless sudo: this is a feature, it makes
   destructive commands deliberate.
7. nginx: use `reload`, never `restart`, and only after `nginx -t`
   passes. If `nginx -t` fails, fix or remove OUR config first.
8. If anything looks wrong mid-deploy, roll back (section 7) rather
   than improvising.

## 1. VM inventory (recon from 2026-09-15, read-only)

Current state, so we know what NOT to touch:

| What                   | Where / port                    | Owner        |
|------------------------|---------------------------------|--------------|
| nginx (reverse proxy)  | 0.0.0.0:80, :443                | system       |
| oligool (FastAPI)      | 127.0.0.1:8001 + :8000 (python) | oligool.service |
| primerool              | 127.0.0.1:8002                  | primerool-serve |
| biochemcalc            | 0.0.0.0:8080                    | biochemcalc-ser |
| PM2 node app           | :3000                           | pm2-kowalski.service |
| svatba RSVP app        | (systemd)                       | svatba.service |
| syncthing              | 127.0.0.1:8384, :22000          | syncthing@kowalski |
| cloudflared tunnels    | 127.0.0.1:20241, :20242        | cloudflared*.service |
| ssh                    | :22                             | system       |

- Ports taken: 22, 53, 80, 443, 3000, 8000, 8001, 8002, 8080, 8384,
  20241, 20242, 22000. Anything else must be verified free before use.
- No docker usable by the deploy user (docker.service runs, but the
  user is not in the docker group and has no sudo). So despite the
  plan mentioning containers, the actual deployment is: native Rust
  static binary + systemd + existing nginx. This is simpler anyway.
- Disk: ~13G free on /. Check before every deploy; each backup is a
  full copy of the data dir.

## 2. bactiment footprint (all of it, nothing outside)

```
/home/kowalski/bactiment/
  app/            # current release binary + static frontend
  versions/       # <timestamp>/ kept releases for rollback
  data/           # SQLite DB + project files (USER DATA)
  backups/        # automatic backups before every push
  logs/           # app logs (journald is primary, files secondary)
```

- Allocated port: 127.0.0.1:8010 (verified free at doc time; re-verify).
- systemd unit: `bactiment.service` (system unit; creation needs sudo,
  run that manually per rule 6).
- nginx: one NEW file `/etc/nginx/sites-available/bactiment` + symlink
  into sites-enabled. We never touch other sites' files.
- Public URL: subdomain or path routed by nginx to 127.0.0.1:8010
  (exact hostname decided at first deploy; cloudflared is an
  alternative if DNS is managed there, but default is nginx).

## 3. Pre-deployment check (every deploy, read-only, safe)

Run all of these and actually read the output:

```bash
ssh proxmox1 '
  echo "=== uptime/load ==="; uptime
  echo "=== disk ==="; df -h / | tail -1
  echo "=== memory ==="; free -h
  echo "=== bactiment port 8010 (must be empty) ==="
  ss -tlnp | grep -E ":8010\b" || echo "free"
  echo "=== foreign services still running ==="
  systemctl is-active oligool svatba pm2-kowalski nginx
  echo "=== nothing new on our other ports ==="
  ss -tlnp
'
```

Abort the deploy if:

- load is high or the box is out of disk/memory,
- something unexpected now listens on 8010 or a port we planned to use,
- any foreign service is down (deploying during an unrelated incident
  only adds confusion),
- `ss` shows new listeners we did not know about (re-plan the port).

## 4. First installation (one time)

On the workstation (this repo):

```bash
cargo build --release
cd frontend && npm ci && npm run build && cd ..
```

Upload to a fresh version dir, never in place:

```bash
ssh proxmox1 'mkdir -p ~/bactiment/versions/$(date +%Y%m%d-%H%M%S) \
                          ~/bactiment/{app,data,backups,logs}'
rsync -av target/release/bactiment-api proxmox1:~/bactiment/versions/<V>/
rsync -av frontend/dist proxmox1:~/bactiment/versions/<V>/
```

Point `app` at the release (a symlink, so rollback is a one-liner):

```bash
ssh proxmox1 'ln -sfn ~/bactiment/versions/<V>/* ~/bactiment/app/'
```

systemd unit `/etc/nginx/...` first the service,
`/etc/systemd/system/bactiment.service` (create with sudo, manually):

```ini
[Unit]
Description=bactiment app
After=network.target

[Service]
User=kowalski
WorkingDirectory=/home/kowalski/bactiment
Environment=BACTIMENT_DATA_DIR=/home/kowalski/bactiment/data
Environment=BACTIMENT_BIND=127.0.0.1:8010
ExecStart=/home/kowalski/bactiment/app/bactiment-api
Restart=on-failure
RestartSec=5

[Install]
WantedBy=multi-user.target
```

```bash
sudo systemctl daemon-reload
sudo systemctl enable --now bactiment.service
```

nginx site `/etc/nginx/sites-available/bactiment` (new file, sudo):

```nginx
server {
    listen 443 ssl;
    server_name bactiment.<domain>;   # decide at deploy time
    # TLS: reuse the certificate strategy of the existing sites
    location / {
        proxy_pass http://127.0.0.1:8010;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        client_max_body_size 512m;    # genome uploads
    }
}
```

```bash
sudo ln -s /etc/nginx/sites-available/bactiment /etc/nginx/sites-enabled/
sudo nginx -t            # MUST pass; if not, remove our file first
sudo systemctl reload nginx   # reload, never restart
```

## 5. Backup policy (before EVERY push of a new version)

The backup always happens before the new binary is activated, and it
must succeed or the deploy aborts:

```bash
ssh proxmox1 '
  set -e
  cd ~/bactiment
  V=$(date +%Y%m%d-%H%M%S)
  # online-safe SQLite backup, then full data dir copy
  sqlite3 data/bactiment.db ".backup backups/db-$V.sqlite" 2>/dev/null \
    || cp data/bactiment.db "backups/db-$V.sqlite"
  tar -C data -czf "backups/data-$V.tar.gz" .
  ls -lh backups/ | tail -3
  df -h / | tail -1
'
```

- Retention: keep the last 10 backups, plus a pinned "known good"
  copy after major data migrations. Older ones are deleted by full
  path inside `~/bactiment/backups` ONLY.
- Backups also run from the app itself nightly (roadmap), this is the
  deploy-time manual guarantee.
- Verify the newest backup is non-trivial in size before continuing.

## 6. Routine update (every later push)

1. Pre-deployment check (section 3).
2. Backup (section 5). No backup, no deploy.
3. Build locally, upload as a NEW version dir, flip the symlink:

```bash
rsync -av target/release/bactiment-api proxmox1:~/bactiment/versions/<V>/
rsync -av frontend/dist proxmox1:~/bactiment/versions/<V>/
ssh proxmox1 'ln -sfn ~/bactiment/versions/<V>/bactiment-api ~/bactiment/app/bactiment-api
              ln -sfn ~/bactiment/versions/<V>/dist ~/bactiment/app/dist'
```

4. Restart ONLY our unit:

```bash
sudo systemctl restart bactiment.service
systemctl status bactiment.service --no-pager
curl -s http://127.0.0.1:8010/api/health
```

5. Smoke test from the browser: create a throwaway project, upload a
  small FASTA, run it, open table + genome views, export a TSV, delete
  the project.
6. Watch logs for a few minutes:

```bash
journalctl -u bactiment.service -f --since "5 min ago"
```

Database schema changes: the app runs migrations at startup and keeps
them additive/backward compatible where possible; the section 5 backup
is the safety net. If a migration is destructive it must be called out
in the release notes and the backup verified extra carefully.

## 7. Rollback

```bash
ssh proxmox1 'ln -sfn ~/bactiment/versions/<PREVIOUS>/bactiment-api ~/bactiment/app/bactiment-api
              ln -sfn ~/bactiment/versions/<PREVIOUS>/dist ~/bactiment/app/dist'
sudo systemctl restart bactiment.service
```

If user data must also be restored (last resort, destroys projects
created after the backup):

```bash
ssh proxmox1 'cd ~/bactiment
              rm -rf data/*          # ONLY inside ~/bactiment/data
              tar -xzf backups/data-<V>.tar.gz -C data'
sudo systemctl restart bactiment.service
```

## 8. Explicitly forbidden (on this VM)

- `kill`, `pkill`, `systemctl restart/stop` for anything not named
  `bactiment*`.
- `docker` anything (not our deployment path, avoid surprises).
- Editing files under other apps' dirs, `/etc/nginx/sites-*` files not
  named `bactiment`, pm2, syncthing folders.
- `rm -rf` with a path that is not strictly inside
  `/home/kowalski/bactiment`.
- Binding ports other than 8010 without re-running the port check.
- `systemctl restart nginx` (reload only, after `nginx -t`).
- Disabling unattended-upgrades or other system services.

## 9. Troubleshooting

| Symptom                          | First move                                   |
|----------------------------------|----------------------------------------------|
| 8010 already bound               | It is either an old bactiment (restart it) or a foreign process (pick another port, update this doc) |
| App up, nginx 502                | Check `journalctl -u bactiment`, then our site file; do not touch other sites |
| Deploy ran, DB locked errors     | Restart only bactiment; SQLite WAL handles the rest |
| Disk full during backup          | Delete oldest backups inside `~/bactiment/backups` only, then re-evaluate |
| Unknown breakage after deploy    | Roll back (section 7), then debug calmly      |

## 10. Open points for first deploy

- Decide the public hostname (nginx server_name) and TLS cert source;
  mirror whatever the existing sites do.
- Confirm MUMmer + BLAST+ versions to install in `~/bactiment/tools`
  (user local, no sudo needed) and pin them.
- Decide backup schedule for the nightly app-side backup (systemd
  timer `bactiment-backup.timer`, sudo, manual).
