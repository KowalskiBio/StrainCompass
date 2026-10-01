#!/usr/bin/env bash
# Deploy the committed, pushed HEAD of this repo to proxmox1.
#
# The one sanctioned way to update production (see deployment.md section 6).
# Every step follows the golden rules: read-only pre-flight first, a
# verified backup before anything changes, the previous version kept for
# rollback, only our own user service touched, and an automatic rollback
# when the new version does not answer its health check.
#
# Usage: deploy/deploy.sh            (from anywhere inside the repo)
set -euo pipefail

HOST=proxmox1
APP='~/straincompass'
SRC='~/straincompass-src'
MIN_FREE_MB=1500
REPO=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
cd "$REPO"

say() { printf '\n=== %s\n' "$*"; }
die() { printf '\nDEPLOY ABORTED: %s\n' "$*" >&2; exit 1; }

# --- 0. deploy only what is committed and pushed -------------------------
say "Checking the working tree"
[ -z "$(git status --porcelain --untracked-files=no)" ] \
  || die "uncommitted changes; commit them first so the server runs a known commit"
git fetch -q origin
[ "$(git rev-parse HEAD)" = "$(git rev-parse '@{u}')" ] \
  || die "HEAD is not pushed (or is behind origin); push first"
COMMIT=$(git rev-parse --short HEAD)
echo "deploying $COMMIT"

# --- 1. pre-flight (read-only) ------------------------------------------
say "Pre-flight on $HOST"
ssh -o ConnectTimeout=15 "$HOST" MIN_FREE_MB=$MIN_FREE_MB 'bash -s' <<'EOF' || die "pre-flight failed"
set -e
load=$(cut -d" " -f1 /proc/loadavg)
free_mb=$(df -Pm / | awk "NR==2{print \$4}")
echo "load $load, free ${free_mb} MB"
awk -v l="$load" "BEGIN{exit !(l < 4)}" || { echo "load too high"; exit 1; }
[ "$free_mb" -ge "$MIN_FREE_MB" ] || { echo "less than $MIN_FREE_MB MB free on /"; exit 1; }
for s in oligool svatba pm2-kowalski nginx; do
  [ "$(systemctl is-active $s)" = active ] || { echo "foreign service $s is not active; not deploying into an incident"; exit 1; }
done
# a restart kills analyses in flight: never deploy under a running one
busy=$(sqlite3 "$HOME/straincompass/data/straincompass.db" "SELECT COUNT(*) FROM runs WHERE status IN ('queued','running')" 2>/dev/null || echo "?")
if [ "$busy" != "0" ]; then
  echo "$busy analysis run(s) queued or running; the restart would kill them. Deploy once they have finished."
  exit 1
fi
owner=$(ss -tlnp 2>/dev/null | awk "/:8010 /" || true)
if [ -n "$owner" ] && ! grep -q straincompass <<<"$owner"; then
  echo "port 8010 is held by something else: $owner"; exit 1
fi
echo "pre-flight OK"
EOF

# --- 2. backup (must succeed) -------------------------------------------
say "Backing up the data"
ssh -o ConnectTimeout=15 "$HOST" "cd $APP && bash -s" <<'EOF' || die "backup failed; nothing was changed"
set -e
V=$(date +%Y%m%d-%H%M%S)
# an incomplete backup is worse than none: it costs disk and looks valid
trap 'st=$?; if [ $st -ne 0 ]; then rm -f "backups/data-$V.tar.gz" "backups/db-$V.sqlite"; echo "removed the incomplete backup $V"; fi' EXIT
sqlite3 data/straincompass.db ".backup backups/db-$V.sqlite" 2>/dev/null \
  || cp data/straincompass.db "backups/db-$V.sqlite"
# a run deleted or written while tar reads it makes tar exit 1; try again
ok=0
for try in 1 2 3; do
  if tar -C data -czf "backups/data-$V.tar.gz" . 2>/dev/null; then ok=1; break; fi
  sleep 5
done
[ $ok = 1 ] || { echo "data kept changing during the backup"; exit 1; }
gzip -t "backups/data-$V.tar.gz"
n_data=$(find data -type f | wc -l)
n_arch=$(tar -tzf "backups/data-$V.tar.gz" | grep -vc '/$')
[ "$n_arch" -ge "$n_data" ] || { echo "archive has $n_arch files, data has $n_data"; exit 1; }
ls -lh "backups/data-$V.tar.gz" "backups/db-$V.sqlite"
EOF

# --- 3. build ------------------------------------------------------------
say "Building the frontend locally"
(cd frontend && npm run build >/dev/null) || die "frontend build failed; nothing was changed"

say "Building the server on $HOST"
rsync -az --delete --exclude target --exclude node_modules --exclude .git \
  --exclude frontend/dist ./ "$HOST:straincompass-src/"
ssh -o ConnectTimeout=15 "$HOST" "cd $SRC && ~/.cargo/bin/cargo build --release 2>&1 | tail -2" \
  || die "server build failed; nothing was changed"
rsync -az --delete frontend/dist/ "$HOST:straincompass/frontend/dist.new/"

# --- 4. swap, health check, automatic rollback ---------------------------
say "Switching to $COMMIT"
ssh -o ConnectTimeout=15 "$HOST" "cd $APP && COMMIT=$COMMIT SRC=$SRC bash -s" <<'EOF'
set -e
P=versions/prev-$(date +%Y%m%d-%H%M%S)
mkdir -p "$P"
cp bin/straincompass-api "$P/"
cp -a frontend/dist "$P/dist"
echo "previous version kept in ~/straincompass/$P"

healthy() {
  for i in $(seq 1 15); do
    sleep 2
    [ "$(curl -s -o /dev/null -w '%{http_code}' http://127.0.0.1:8010/api/health)" = 200 ] && return 0
  done
  return 1
}

# everything that can fail slowly happens before the stop; while the
# service is down only renames run, and any error still restarts it
cp "${SRC/#\~/$HOME}/target/release/straincompass-api" bin/straincompass-api.new
rm -rf frontend/dist.old
systemctl --user stop straincompass
trap 'systemctl --user start straincompass; echo "error during the switch; service restarted" >&2' ERR
mv bin/straincompass-api.new bin/straincompass-api
mv frontend/dist frontend/dist.old
mv frontend/dist.new frontend/dist
systemctl --user start straincompass
trap - ERR
if healthy; then
  echo "$COMMIT $(date -Is)" >> deployed.log
  echo "healthy: $COMMIT is live"
  exit 0
fi

echo "NEW VERSION FAILED ITS HEALTH CHECK - rolling back to $P" >&2
systemctl --user stop straincompass || true
cp "$P/straincompass-api" bin/straincompass-api
rm -rf frontend/dist
cp -a "$P/dist" frontend/dist
systemctl --user start straincompass
healthy && echo "rolled back, the previous version is running" >&2 \
  || echo "ROLLBACK ALSO UNHEALTHY - check: journalctl --user -u straincompass" >&2
exit 1
EOF

say "Done: $COMMIT is live on https://sc.ubch.sci.muni.cz"
