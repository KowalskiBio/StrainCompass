import { useCallback, useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";
import { api } from "../api";
import type { Project } from "../types";
import { Button, EmptyState, Spinner } from "../components/ui";
import { NewProjectDialog } from "./NewProjectDialog";

export default function ProjectsPage() {
  const navigate = useNavigate();
  const [projects, setProjects] = useState<Project[]>([]);
  const [loading, setLoading] = useState(true);
  const [newOpen, setNewOpen] = useState(false);
  const [search, setSearch] = useState("");
  const [sortDesc, setSortDesc] = useState(true);

  const reload = useCallback(() => {
    api
      .listProjects()
      .then(setProjects)
      .catch(() => {})
      .finally(() => setLoading(false));
  }, []);

  useEffect(reload, [reload]);

  const filtered = projects
    .filter(
      (p) =>
        !search ||
        p.name.toLowerCase().includes(search.toLowerCase()) ||
        p.organism.toLowerCase().includes(search.toLowerCase()),
    )
    .sort((a, b) =>
      sortDesc
        ? b.created_at.localeCompare(a.created_at)
        : a.created_at.localeCompare(b.created_at),
    );

  return (
    <div className="max-w-5xl mx-auto px-4 py-8">
      <div className="flex flex-wrap items-center justify-between gap-4 mb-6">
        <div>
          <h1 className="text-2xl font-semibold tracking-tight">Projects</h1>
          <p className="text-zinc-500 mt-1">
            Each project holds one annotated reference genome and the query
            genomes you compare against it.
          </p>
        </div>
        <button
          className="h-11 px-5 rounded-lg bg-zinc-900 text-white text-[15px] hover:bg-zinc-700"
          onClick={() => setNewOpen(true)}
        >
          New project
        </button>
      </div>

      <div className="flex items-center justify-between mb-4 gap-4">
        <input
          type="search"
          placeholder="Search projects"
          value={search}
          onChange={(e) => setSearch(e.target.value)}
          className="h-11 w-72 px-3 rounded-lg border border-zinc-300 bg-white text-[15px]"
        />
        <button
          className="text-sm text-zinc-500 hover:text-zinc-900 whitespace-nowrap"
          onClick={() => setSortDesc(!sortDesc)}
        >
          Sorted by date {sortDesc ? "(newest first)" : "(oldest first)"}
        </button>
      </div>

      {loading && projects.length === 0 ? (
        <div className="flex items-center gap-3 text-zinc-500 py-16 justify-center">
          <Spinner /> Loading projects...
        </div>
      ) : filtered.length === 0 ? (
        <EmptyState
          title={search ? "No projects match your search" : "No projects yet"}
          hint={
            search
              ? "Try a different search term, or clear it to see all projects."
              : "Create your first project: upload a reference genome with its annotation, add query genomes, and run the comparison."
          }
          action={
            !search && (
              <Button size="lg" onClick={() => setNewOpen(true)}>
                Create a project
              </Button>
            )
          }
        />
      ) : (
        <ul className="divide-y divide-zinc-200 border border-zinc-200 rounded-xl bg-white">
          {filtered.map((p) => (
            <li
              key={p.id}
              className="flex items-center justify-between gap-4 px-5 py-4 hover:bg-zinc-50 cursor-pointer"
              onClick={() => navigate(`/projects/${p.id}`)}
            >
              <div className="min-w-0">
                <p className="font-medium truncate">{p.name}</p>
                <p className="text-sm text-zinc-500 truncate">
                  {p.organism && p.organism !== "bacteria" ? `${p.organism} - ` : ""}
                  {p.has_reference ? "reference ready" : "no reference yet"}
                </p>
              </div>
              <div className="flex items-center gap-6 shrink-0 text-sm text-zinc-500">
                <span>
                  {p.n_queries} quer{p.n_queries === 1 ? "y" : "ies"}
                </span>
                <span>
                  {p.n_runs} run{p.n_runs === 1 ? "" : "s"}
                </span>
                <span className="hidden sm:inline">
                  created {formatRelative(p.created_at)}
                </span>
              </div>
            </li>
          ))}
        </ul>
      )}

      {newOpen && (
        <NewProjectDialog
          onClose={() => setNewOpen(false)}
          onCreated={(id) => {
            setNewOpen(false);
            navigate(`/projects/${id}`);
          }}
        />
      )}
    </div>
  );
}

function formatRelative(iso: string): string {
  const t = new Date(iso + (iso.endsWith("Z") ? "" : "Z")).getTime();
  if (isNaN(t)) return iso;
  const min = Math.round((Date.now() - t) / 60000);
  if (min < 1) return "just now";
  if (min < 60) return `${min} minute${min > 1 ? "s" : ""} ago`;
  const h = Math.round(min / 60);
  if (h < 24) return `${h} hour${h > 1 ? "s" : ""} ago`;
  const d = Math.round(h / 24);
  if (d < 30) return `${d} day${d > 1 ? "s" : ""} ago`;
  return new Date(iso).toLocaleDateString();
}
