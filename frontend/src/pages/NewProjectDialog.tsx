import { useState } from "react";
import { api } from "../api";
import { Button, ErrorBox, Modal } from "../components/ui";

/** Small dialog that creates a project (name + optional organism). */
export function NewProjectDialog({
  onClose,
  onCreated,
}: {
  onClose: () => void;
  onCreated: (id: number) => void;
}) {
  const [name, setName] = useState("");
  const [organism, setOrganism] = useState("");
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);

  async function create() {
    if (!name.trim()) {
      setError("Please give the project a name.");
      return;
    }
    setBusy(true);
    setError(null);
    try {
      const p = await api.createProject(name.trim(), organism.trim() || undefined);
      onCreated(p.id);
    } catch (e) {
      setError((e as Error).message);
      setBusy(false);
    }
  }

  return (
    <Modal
      open
      onClose={onClose}
      title="New project"
    >
      <div className="space-y-4">
        {error && <ErrorBox message={error} />}
        <label className="block">
          <span className="text-sm font-medium text-zinc-700">Project name</span>
          <input
            autoFocus
            value={name}
            onChange={(e) => setName(e.target.value)}
            onKeyDown={(e) => e.key === "Enter" && create()}
            placeholder="e.g. Listeria monocytogenes comparison"
            className="mt-1.5 w-full h-11 px-3 rounded-lg border border-zinc-300 text-[15px] focus:border-zinc-500 outline-none"
          />
        </label>
        <label className="block">
          <span className="text-sm font-medium text-zinc-700">
            Organism <span className="text-zinc-400 font-normal">(optional)</span>
          </span>
          <input
            value={organism}
            onChange={(e) => setOrganism(e.target.value)}
            onKeyDown={(e) => e.key === "Enter" && create()}
            placeholder="e.g. L. monocytogenes EGD-e"
            className="mt-1.5 w-full h-11 px-3 rounded-lg border border-zinc-300 text-[15px] focus:border-zinc-500 outline-none"
          />
        </label>
        <p className="text-sm text-zinc-500">
          After creating the project you will add the reference genome and the
          query genomes.
        </p>
        <div className="flex justify-end gap-3 pt-2">
          <Button variant="ghost" onClick={onClose}>
            Cancel
          </Button>
          <Button onClick={create} disabled={busy}>
            {busy ? "Creating..." : "Create project"}
          </Button>
        </div>
      </div>
    </Modal>
  );
}
