# StrainCompass deployment status (2026-09-15)

The app is built, deployed and running. One manual step remains that needs
your sudo password.

## Running now

- Service: `systemctl --user status straincompass` on proxmox1
- Listens on: `0.0.0.0:8010` (deliberately public, see deployment.md
  section 11 - reachable at http://147.251.22.180:8010/ )
- App home: `~/straincompass/` (bin, frontend/dist, data)
- Tools (user-local on the VM): `~/.local/opt/MUMmer3.23`, `~/.local/opt/blast`
- Source used for the build: `~/straincompass-src` (matches git commit e7ebfa5 + this commit)
- Verified end to end on the VM: project, uploads, run with real
  nucmer/dnadiff/blastn, all result endpoints, exports, SPA serving.

## Access right now (no sudo, no DNS)

    http://147.251.22.180:8010/

Works from the campus network and over the VPN. If a colleague outside
the university cannot reach it, the perimeter firewall is filtering the
high port - use the port 80 route in deployment.md section 11, which
needs one sudo step and then serves http://147.251.22.180/ .

An SSH tunnel still works too, if you want it private again for a moment:

    ssh -L 8010:127.0.0.1:8010 proxmox1   # then http://localhost:8010

Note: the app has no authentication. Anyone who can reach the port can
create projects, upload 600 MB files and start runs on this VM.

## The one sudo step (public URL + TLS)

The DNS record now exists. The granted hostname is **`sc.ubch.sci.muni.cz`**
— that is the name the request actually asked for, not
`straincompass.ubch.sci.muni.cz`, which this file previously assumed and which
was never created. It is a CNAME to `oligool.biochem.muni.cz` ->
147.251.22.180, created 2026-09-15, same shape as the sibling sites.

Verified: `http://sc.ubch.sci.muni.cz/` already reaches the VM on port 80
and returns nginx's 404, because no vhost claims that Host header yet.
Port 80 is open through the perimeter, so the Let's Encrypt HTTP-01
challenge will work.

The site file is prepared at `deploy/nginx-straincompass.conf` and staged on the
VM at `~/straincompass/nginx-sc.conf`. The deploy user has no passwordless sudo
(golden rule 6), so these three run manually:

    sudo cp ~/straincompass/nginx-sc.conf /etc/nginx/sites-available/sc.ubch.sci.muni.cz
    sudo ln -s /etc/nginx/sites-available/sc.ubch.sci.muni.cz /etc/nginx/sites-enabled/
    sudo nginx -t && sudo systemctl reload nginx

Then issue the certificate (certbot 2.9.0 is already installed, and it
rewrites the file in place to add the TLS block and the 80 -> 443 redirect):

    sudo certbot --nginx -d sc.ubch.sci.muni.cz

Then open https://sc.ubch.sci.muni.cz

To undo the whole thing:

    sudo rm /etc/nginx/sites-enabled/sc.ubch.sci.muni.cz
    sudo nginx -t && sudo systemctl reload nginx

## After TLS works: close the direct port

Until this is done the app is still reachable, unencrypted and
unauthenticated, at `http://147.251.22.180:8010/`. Once HTTPS is confirmed:

    sed -i 's|STRAINCOMPASS_BIND=0.0.0.0:8010|STRAINCOMPASS_BIND=127.0.0.1:8010|' \
        ~/.config/systemd/user/straincompass.service
    systemctl --user daemon-reload && systemctl --user restart straincompass

No sudo needed (it is a user unit). This also retires golden rule 3's
exception and deployment.md section 11.

**TLS is not authentication.** The API still has none: anyone who reaches
https://sc.ubch.sci.muni.cz can create and delete projects, upload 600 MB
per request and start nucmer/dnadiff/blastn runs on this shared production
VM. HTTPS stops network eavesdropping, nothing else. A public hostname
makes the app easier to find than a bare IP on a high port, so this is now
the top open item.

## Updating later

    # on the workstation
    rsync -az --delete --exclude target --exclude node_modules --exclude .git \
        /home/kowalski/Work/straincompass/ proxmox1:straincompass-src/
    # on the VM
    cd ~/straincompass-src && ~/.cargo/bin/cargo build --release
    systemctl --user stop straincompass
    cp target/release/straincompass-api ~/straincompass/bin/
    # frontend (if changed): rebuild locally and
    rsync -az frontend/dist/ proxmox1:straincompass/frontend/dist/
    systemctl --user start straincompass

Back up `~/straincompass/data` before deploys that change the database schema
(currently: the whole app state lives there, SQLite + project files).
