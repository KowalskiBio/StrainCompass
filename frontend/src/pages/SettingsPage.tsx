import { useEffect, useState } from "react";
import { api } from "../api";
import { Button, ErrorBox } from "../components/ui";

export default function SettingsPage() {
  const [hasKey, setHasKey] = useState(false);
  const [masked, setMasked] = useState<string | null>(null);
  const [key, setKey] = useState("");
  const [error, setError] = useState<string | null>(null);
  const [notice, setNotice] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const [loaded, setLoaded] = useState(false);

  useEffect(() => {
    api
      .getSettings()
      .then((s) => {
        setHasKey(s.has_ncbi_api_key);
        setMasked(s.ncbi_api_key);
      })
      .catch((e) => setError((e as Error).message))
      .finally(() => setLoaded(true));
  }, []);

  async function save() {
    const trimmed = key.trim();
    if (!trimmed) {
      setError("Paste the API key first.");
      return;
    }
    setBusy(true);
    setError(null);
    setNotice(null);
    try {
      const r = await api.putNcbiKey(trimmed);
      setHasKey(true);
      setMasked(r.masked);
      setKey("");
      setNotice("The API key has been saved.");
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy(false);
    }
  }

  async function remove() {
    setBusy(true);
    setError(null);
    setNotice(null);
    try {
      await api.deleteNcbiKey();
      setHasKey(false);
      setMasked(null);
      setNotice("The API key has been removed.");
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy(false);
    }
  }

  return (
    <div className="max-w-2xl mx-auto px-4 py-8">
      <h1 className="text-2xl font-semibold tracking-tight">Settings</h1>

      <section className="mt-6 border border-zinc-200 rounded-xl bg-white p-6 dark:border-zinc-800 dark:bg-zinc-900">
        <h2 className="text-lg font-medium">NCBI API key (optional)</h2>
        <p className="text-[15px] text-zinc-500 mt-1 dark:text-zinc-400">
          An API key raises the limit on NCBI downloads, which speeds up
          fetching reference genomes by accession. Without a key everything
          still works, just slower.
        </p>

        {error && (
          <div className="mt-4">
            <ErrorBox message={error} />
          </div>
        )}
        {notice && (
          <p className="mt-4 text-[15px] text-emerald-700 dark:text-emerald-400">{notice}</p>
        )}

        {loaded && hasKey && (
          <div className="mt-4 flex items-center justify-between gap-4 rounded-lg border border-emerald-200 bg-emerald-50 px-4 py-3 dark:border-emerald-900 dark:bg-emerald-950/40">
            <p className="text-[15px] text-emerald-900 dark:text-emerald-300">
              A key is saved on the server: <span className="font-mono">{masked}</span>
            </p>
            <Button variant="danger" onClick={remove} disabled={busy}>
              Remove
            </Button>
          </div>
        )}

        <div className="mt-4">
          <label className="block">
            <span className="text-sm font-medium text-zinc-700 dark:text-zinc-300">
              {hasKey ? "Replace the saved key" : "Add a key"}
            </span>
            <div className="mt-1.5 flex gap-2">
              <input
                value={key}
                onChange={(e) => setKey(e.target.value)}
                placeholder="e.g. a1b2c3d4e5f6g7h8i9j0"
                className="flex-1 h-11 px-3 rounded-lg border border-zinc-300 text-[15px] font-mono focus:border-zinc-500 outline-none dark:border-zinc-700 dark:bg-zinc-900"
              />
              <Button onClick={save} disabled={busy}>
                {busy ? "Saving..." : "Save"}
              </Button>
            </div>
          </label>
          <p className="text-xs text-zinc-400 mt-2 dark:text-zinc-500">
            Get a key at the NCBI website (Account, API Key Management). The
            key is stored on this server only.
          </p>
        </div>
      </section>
    </div>
  );
}
