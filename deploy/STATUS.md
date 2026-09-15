# Bactiment deployment status (2026-09-15)

The app is built, deployed and running. One manual step remains that needs
your sudo password.

## Running now

- Service: `systemctl --user status bactiment` on proxmox1
- Listens on: `127.0.0.1:8010` (localhost only, as required)
- App home: `~/bactiment/` (bin, frontend/dist, data)
- Tools (user-local on the VM): `~/.local/opt/MUMmer3.23`, `~/.local/opt/blast`
- Source used for the build: `~/bactiment-src` (matches git commit e7ebfa5 + this commit)
- Verified end to end on the VM: project, uploads, run with real
  nucmer/dnadiff/blastn, all result endpoints, exports, SPA serving.

## Immediate access without any sudo

From your workstation:

    ssh -L 8010:127.0.0.1:8010 proxmox1

then open http://localhost:8010

## The one sudo step (public URL)

The VM nginx serves one subdomain per app (biocal, primerool, oligool).
Bactiment expects `bactiment.ubch.sci.muni.cz`. The site file is prepared at
`deploy/nginx-bactiment.conf` (also copied to the VM at
`~/bactiment/nginx-bactiment.conf`). Steps:

1. Create the DNS record `bactiment.ubch.sci.muni.cz` pointing at the VM
   (wherever your DNS is managed).
2. On the VM:

       sudo cp ~/bactiment/nginx-bactiment.conf /etc/nginx/sites-available/bactiment.ubch.sci.muni.cz
       sudo ln -s /etc/nginx/sites-available/bactiment.ubch.sci.muni.cz /etc/nginx/sites-enabled/
       sudo nginx -t && sudo systemctl reload nginx
       sudo certbot --nginx -d bactiment.ubch.sci.muni.cz

3. Open https://bactiment.ubch.sci.muni.cz

The conf file proxies everything to 127.0.0.1:8010 with a 600 MB body limit
(same limit as the app's upload limit). Certbot rewrites the file to add the
TLS server block, same as the other sites.

If you prefer a different subdomain, edit `server_name` in the conf file.

## Updating later

    # on the workstation
    rsync -az --delete --exclude target --exclude node_modules --exclude .git \
        /home/kowalski/Work/bactiment/ proxmox1:bactiment-src/
    # on the VM
    cd ~/bactiment-src && ~/.cargo/bin/cargo build --release
    systemctl --user stop bactiment
    cp target/release/bactiment-api ~/bactiment/bin/
    # frontend (if changed): rebuild locally and
    rsync -az frontend/dist/ proxmox1:bactiment/frontend/dist/
    systemctl --user start bactiment

Back up `~/bactiment/data` before deploys that change the database schema
(currently: the whole app state lives there, SQLite + project files).
