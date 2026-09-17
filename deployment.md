# straincompass deployment guide (VM: proxmox1)

This VM is PRODUCTION and hosts other people's apps. Every deployment
follows this document. If any step would conflict with another app,
STOP and reconsider instead of forcing it.

## 0. Golden safety rules (non negotiable)

1. NEVER stop, restart, kill or "fix" any process that does not belong
   to straincompass. No `kill` on foreign PIDs, no `systemctl restart` of
   foreign units, no `pm2 restart`, nothing.
2. NEVER edit, overwrite or delete another app's files, nginx site
   configs, or data directories. straincompass gets its own files only.
3. NEVER bind a port that is already in use. straincompass currently binds
   0.0.0.0:8010 on purpose (see section 11): the owner asked for direct
   http://IP:port access while the DNS record is pending. Do NOT "fix"
   this back to 127.0.0.1 without asking. Every OTHER port stays
   localhost-only, and the preferred public side is still nginx.
4. ALWAYS back up straincompass's data before every deploy/push of a new
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

## 2. straincompass footprint (all of it, nothing outside)

```
/home/kowalski/straincompass/
  app/            # current release binary + static frontend
  versions/       # <timestamp>/ kept releases for rollback
  data/           # SQLite DB + project files (USER DATA)
  backups/        # automatic backups before every push
  logs/           # app logs (journald is primary, files secondary)
```

- Allocated port: 8010, currently bound on 0.0.0.0 (section 11).
- systemd unit: `straincompass.service` (system unit; creation needs sudo,
  run that manually per rule 6).
- nginx: one NEW file `/etc/nginx/sites-available/straincompass` + symlink
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
  echo "=== straincompass port 8010 (must be empty) ==="
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
ssh proxmox1 'mkdir -p ~/straincompass/versions/$(date +%Y%m%d-%H%M%S) \
                          ~/straincompass/{app,data,backups,logs}'
rsync -av target/release/straincompass-api proxmox1:~/straincompass/versions/<V>/
rsync -av frontend/dist proxmox1:~/straincompass/versions/<V>/
```

Point `app` at the release (a symlink, so rollback is a one-liner):

```bash
ssh proxmox1 'ln -sfn ~/straincompass/versions/<V>/* ~/straincompass/app/'
```

systemd unit `/etc/nginx/...` first the service,
`/etc/systemd/system/straincompass.service` (create with sudo, manually):

```ini
[Unit]
Description=straincompass app
After=network.target

[Service]
User=kowalski
WorkingDirectory=/home/kowalski/straincompass
Environment=STRAINCOMPASS_DATA_DIR=/home/kowalski/straincompass/data
Environment=STRAINCOMPASS_BIND=0.0.0.0:8010   # see section 11; 127.0.0.1:8010 once DNS + nginx are in place
ExecStart=/home/kowalski/straincompass/app/straincompass-api
Restart=on-failure
RestartSec=5

[Install]
WantedBy=multi-user.target
```

```bash
sudo systemctl daemon-reload
sudo systemctl enable --now straincompass.service
```

nginx site `/etc/nginx/sites-available/straincompass` (new file, sudo):

```nginx
server {
    listen 443 ssl;
    server_name straincompass.<domain>;   # decide at deploy time
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
sudo ln -s /etc/nginx/sites-available/straincompass /etc/nginx/sites-enabled/
sudo nginx -t            # MUST pass; if not, remove our file first
sudo systemctl reload nginx   # reload, never restart
```

## 5. Backup policy (before EVERY push of a new version)

The backup always happens before the new binary is activated, and it
must succeed or the deploy aborts:

```bash
ssh proxmox1 '
  set -e
  cd ~/straincompass
  V=$(date +%Y%m%d-%H%M%S)
  # online-safe SQLite backup, then full data dir copy
  sqlite3 data/straincompass.db ".backup backups/db-$V.sqlite" 2>/dev/null \
    || cp data/straincompass.db "backups/db-$V.sqlite"
  tar -C data -czf "backups/data-$V.tar.gz" .
  ls -lh backups/ | tail -3
  df -h / | tail -1
'
```

- Retention: keep the last 10 backups, plus a pinned "known good"
  copy after major data migrations. Older ones are deleted by full
  path inside `~/straincompass/backups` ONLY.
- Backups also run from the app itself nightly (roadmap), this is the
  deploy-time manual guarantee.
- Verify the newest backup is non-trivial in size before continuing.

## 6. Routine update (every later push)

1. Pre-deployment check (section 3).
2. Backup (section 5). No backup, no deploy.
3. Build locally, upload as a NEW version dir, flip the symlink:

```bash
rsync -av target/release/straincompass-api proxmox1:~/straincompass/versions/<V>/
rsync -av frontend/dist proxmox1:~/straincompass/versions/<V>/
ssh proxmox1 'ln -sfn ~/straincompass/versions/<V>/straincompass-api ~/straincompass/app/straincompass-api
              ln -sfn ~/straincompass/versions/<V>/dist ~/straincompass/app/dist'
```

4. Restart ONLY our unit:

```bash
sudo systemctl restart straincompass.service
systemctl status straincompass.service --no-pager
curl -s http://127.0.0.1:8010/api/health
```

5. Smoke test from the browser: create a throwaway project, upload a
  small FASTA, run it, open table + genome views, export a TSV, delete
  the project.
6. Watch logs for a few minutes:

```bash
journalctl -u straincompass.service -f --since "5 min ago"
```

Database schema changes: the app runs migrations at startup and keeps
them additive/backward compatible where possible; the section 5 backup
is the safety net. If a migration is destructive it must be called out
in the release notes and the backup verified extra carefully.

## 7. Rollback

```bash
ssh proxmox1 'ln -sfn ~/straincompass/versions/<PREVIOUS>/straincompass-api ~/straincompass/app/straincompass-api
              ln -sfn ~/straincompass/versions/<PREVIOUS>/dist ~/straincompass/app/dist'
sudo systemctl restart straincompass.service
```

If user data must also be restored (last resort, destroys projects
created after the backup):

```bash
ssh proxmox1 'cd ~/straincompass
              rm -rf data/*          # ONLY inside ~/straincompass/data
              tar -xzf backups/data-<V>.tar.gz -C data'
sudo systemctl restart straincompass.service
```

## 8. Explicitly forbidden (on this VM)

- `kill`, `pkill`, `systemctl restart/stop` for anything not named
  `straincompass*`.
- `docker` anything (not our deployment path, avoid surprises).
- Editing files under other apps' dirs, `/etc/nginx/sites-*` files not
  named `straincompass`, pm2, syncthing folders.
- `rm -rf` with a path that is not strictly inside
  `/home/kowalski/straincompass`.
- Binding ports other than 8010 without re-running the port check.
- Exposing any port other than 8010 on 0.0.0.0.
- `systemctl restart nginx` (reload only, after `nginx -t`).
- Disabling unattended-upgrades or other system services.

## 9. Troubleshooting

| Symptom                          | First move                                   |
|----------------------------------|----------------------------------------------|
| 8010 already bound               | It is either an old straincompass (restart it) or a foreign process (pick another port, update this doc) |
| App up, nginx 502                | Check `journalctl -u straincompass`, then our site file; do not touch other sites |
| Deploy ran, DB locked errors     | Restart only straincompass; SQLite WAL handles the rest |
| Disk full during backup          | Delete oldest backups inside `~/straincompass/backups` only, then re-evaluate |
| Unknown breakage after deploy    | Roll back (section 7), then debug calmly      |

## 10. Open points for first deploy

- Decide the public hostname (nginx server_name) and TLS cert source;
  mirror whatever the existing sites do.
- Prodigal is an optional seventh tool, found the same way (PATH or
  `STRAINCOMPASS_TOOLS_DIRS`). It only predicts the genes inside gained
  regions; when it is absent every run still succeeds and the gained
  tables simply report no gene counts. It also adds two small per-query
  artifacts, `gained_regions.tsv` and `gained_regions.fa`, the latter
  bounded by the divergent fraction of the query genome.
- Confirm MUMmer + BLAST+ versions to install in `~/straincompass/tools`
  (user local, no sudo needed) and pin them.
- Decide backup schedule for the nightly app-side backup (systemd
  timer `straincompass-backup.timer`, sudo, manual).

## 11. Temporary public access by IP (until DNS exists)

The owner asked for the app to be reachable from other computers before
DNS existed, so straincompass is deliberately exposed. The hostname
`sc.ubch.sci.muni.cz` resolves as of 2026-09-15, so this section is now
historical: it is the record of that decision, to be removed once the
vhost + certificate are installed and the app is re-bound to localhost
(see deploy/STATUS.md).

Reality check, because it differs from sections 2 and 4:

- The running unit is a **user** unit, `~/.config/systemd/user/straincompass.service`,
  not the system unit sketched in section 4. `systemctl --user ...`,
  and `loginctl show-user kowalski -p Linger` is `Linger=yes`, so it
  survives reboot without anyone logging in.
- It binds `0.0.0.0:8010` (`STRAINCOMPASS_BIND` in that unit).
- The VM's address is `147.251.22.180` (ens18).

Two ways in, in order of preference:

1. **Port 80 via nginx** — `http://147.251.22.180/`. Needs the one sudo
   step below. Preferred because port 80 is already open through the
   university perimeter (the sibling sites work from anywhere), whereas
   a high port like 8010 may be filtered outside the campus network.

   ```bash
   sudo cp ~/straincompass/nginx-straincompass-ip.conf /etc/nginx/sites-available/straincompass-ip
   sudo ln -s /etc/nginx/sites-available/straincompass-ip /etc/nginx/sites-enabled/
   sudo nginx -t && sudo systemctl reload nginx
   ```

   The file matches `server_name 147.251.22.180` only, is not a
   `default_server`, and touches no other site, so biocal / oligool /
   primerool are unaffected. To undo:

   ```bash
   sudo rm /etc/nginx/sites-enabled/straincompass-ip
   sudo nginx -t && sudo systemctl reload nginx
   ```

2. **Direct port** — `http://147.251.22.180:8010/`. Already live, no
   sudo needed. Verified working from the campus/VPN network; whether it
   passes the perimeter firewall from the open internet is unconfirmed.

### What exposure actually means here

The API has **no authentication of any kind**. Anyone who can reach the
port can create and delete projects, upload up to 600 MB per request,
and start nucmer/dnadiff/blastn runs on this shared production VM. The
NCBI API key is masked in `GET /settings`, so that is not leaking, but
compute and disk are open. Acceptable for a short testing window with
colleagues; not acceptable as a permanent state.

### Going back to private, once DNS is ready

```bash
# 1. real vhost + TLS (deploy/nginx-straincompass.conf, section 4)
# 2. drop the IP vhost (undo command above)
# 3. re-bind to localhost
sed -i 's|STRAINCOMPASS_BIND=0.0.0.0:8010|STRAINCOMPASS_BIND=127.0.0.1:8010|' \
    ~/.config/systemd/user/straincompass.service
systemctl --user daemon-reload && systemctl --user restart straincompass
# 4. flip golden rule 3 back to "localhost only"
```

### The "operation is insecure" crash (fixed 2026-09-15)

Symptom: Zen (a Firefox fork) showed the app's error card with
`SecurityError: The operation is insecure.` Chromium was fine.

This was NOT caused by the public IP, by HTTP vs HTTPS, or by the
service being a user unit. It was a frontend render loop:

- `ProjectPage` passed `onStateChange` to `ResultsTables` as an inline
  arrow, so its identity changed on every parent render.
- `ResultsTables` listed that callback in a `useEffect` dependency
  array, so the effect re-ran on every parent render, calling
  `setParam` -> `setSearchParams` -> re-render -> repeat.
- Each `setParam` is one `history.replaceState`. Measured: **202 history
  writes for a single tab switch.**
- Gecko caps history writes at roughly 200 per 10 seconds and throws
  `SecurityError` past the cap. Chromium only throttles silently, which
  is why it never surfaced during development.

Fix: the callback lives in a ref in `ResultsTables`, the handlers in
`ProjectPage` are `useCallback`-stable, and `setParamsBatch` writes
several query params in one history entry (the genome view's
`onRangeChange` wrote three in a row on every pan/zoom). Same tab switch
now costs **3 writes**, and switching table kind costs 1, with the URL
still tracking state so deep links keep working.

Worth remembering: a browser that "works" is not proof of correctness
here. Chromium hid a 200x write amplification that only a stricter
engine reported.

## 12. Disk footprint

`dnadiff` writes `<prefix>.snps`, a full SNP listing of ~19 MB per
bacterial genome pair. Nothing reads it: not the engine, not the API,
and it is not in the Files panel whitelist in
`crates/api/src/routes/runfiles.rs`. On the first real project it was
304.6 MB of a 440 MB footprint - 69% of the project, invisible and
permanent.

`Tools::run_dnadiff` now deletes it as soon as the report exists
(`crates/engine/src/tools.rs`). The same 16-query Listeria project went
from 442 MB to 137 MB on disk. Future runs of that shape cost ~80 MB
instead of ~386 MB.

What a run legitimately keeps per query: `cmp.delta` (parsed for
coverage and gaps), `cmp.report` (shown and downloadable), the result
TSVs, and the BLAST database for the gene panel.

Still open: `spawn_maintenance` in `crates/api/src/main.rs` is an empty
loop with a placeholder comment, so nothing reclaims anything on a
schedule. Old runs and deleted projects are only cleaned when a user
deletes them. Worth revisiting before the VM gets tight - it is shared
with other people's apps, and `/` sits around 75%.

Note the app buffers each uploaded file fully in RAM (`field.bytes()`),
so the 600 MB per-request body limit is also a per-request memory
budget, and the unit has `MemoryMax=infinity`.

## 13. Delete looked broken but always worked (fixed 2026-09-15)

Symptom: deleting a project showed `Unexpected token 'd', "deleted" is
not valid JSON`, the dialog stayed open, and a second click produced
`This project does not exist (anymore).`

Cause: three handlers answered with a bare `"deleted"` string
(`projects.rs`, `runs.rs`, `uploads.rs` delete routes) while the frontend
`request()` helper always called `resp.json()`. The server had already
done the delete; only the response parsing failed. So `onDone` never ran,
the dialog never closed, and the second attempt hit a row that was
already gone.

It affected delete project, delete run AND delete file. Anyone using it
would conclude nothing could be deleted, while everything was in fact
being deleted.

Fixed on both sides: the three routes return `{"status":"deleted"}` with
`content-type: application/json`, and `request()` now only parses when the
server actually says it sent JSON (204 and non-JSON bodies return
undefined instead of throwing).

Lesson for this codebase: a handler whose success body is not JSON is a
silent trap, because the failure surfaces as a client-side parse error
that looks like the operation failed.

## 14. Run names (added 2026-09-15)

Runs can be labelled from the Runs tab, so several comparisons of the same
project stay distinguishable ("strict 95% coverage" vs "with panel
recheck") instead of being "Run #4" and "Run #5".

- Schema: `runs.name TEXT`, nullable. Added by the first additive
  migration in this codebase - `add_column_if_missing` in
  `crates/api/src/db.rs` checks `pragma table_info` on every start,
  because `CREATE TABLE IF NOT EXISTS` does nothing to an existing table
  and SQLite has no `ADD COLUMN IF NOT EXISTS`. Add future columns the
  same way.
- API: `PUT /runs/{id}/name` with `{"name": "..."}`. An empty or
  whitespace-only name clears the label; over 120 characters is a 400.
- UI: `runLabel()` in `frontend/src/types.ts` is the single place that
  decides what a run is called, falling back to `Run #<id>`. A named run
  still shows its `#id` so it stays identifiable in logs and file paths.
