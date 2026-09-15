import { type ReactNode, useEffect, useRef } from "react";
import type { Call } from "../types";

export function Button({
  children,
  variant = "primary",
  size = "md",
  className = "",
  ...props
}: {
  children: ReactNode;
  variant?: "primary" | "secondary" | "ghost" | "danger";
  size?: "md" | "lg";
} & React.ButtonHTMLAttributes<HTMLButtonElement>) {
  const base =
    "inline-flex items-center justify-center gap-2 font-medium rounded-lg transition-colors disabled:opacity-50 disabled:pointer-events-none select-none";
  const sizes = size === "lg" ? "h-12 px-6 text-base" : "h-11 px-4 text-[15px]";
  const variants = {
    primary:
      "bg-zinc-900 text-white hover:bg-zinc-700 dark:bg-zinc-100 dark:text-zinc-900 dark:hover:bg-zinc-300",
    secondary:
      "bg-white text-zinc-900 border border-zinc-300 hover:bg-zinc-100 dark:bg-zinc-900 dark:text-zinc-100 dark:border-zinc-700 dark:hover:bg-zinc-800",
    ghost:
      "text-zinc-600 hover:text-zinc-900 hover:bg-zinc-100 dark:text-zinc-400 dark:hover:text-zinc-100 dark:hover:bg-zinc-800",
    danger:
      "bg-white text-red-700 border border-red-300 hover:bg-red-50 dark:bg-zinc-900 dark:text-red-400 dark:border-red-900 dark:hover:bg-red-950/40",
  }[variant];
  return (
    <button className={`${base} ${sizes} ${variants} ${className}`} {...props}>
      {children}
    </button>
  );
}

export function CallBadge({ call }: { call: Call }) {
  const map = {
    PRESENT:
      "bg-emerald-50 text-emerald-800 border-emerald-200 dark:bg-emerald-950/40 dark:text-emerald-300 dark:border-emerald-900",
    PARTIAL:
      "bg-amber-50 text-amber-800 border-amber-200 dark:bg-amber-950/40 dark:text-amber-300 dark:border-amber-900",
    ABSENT:
      "bg-zinc-100 text-zinc-600 border-zinc-200 dark:bg-zinc-800 dark:text-zinc-400 dark:border-zinc-700",
  } as const;
  const icon = { PRESENT: "\u2713", PARTIAL: "\u25D0", ABSENT: "\u2715" }[call];
  return (
    <span
      className={`inline-flex items-center gap-1.5 px-2 py-0.5 rounded-full border text-xs font-semibold whitespace-nowrap ${map[call]}`}
    >
      <span aria-hidden>{icon}</span>
      {call === "PRESENT" ? "Present" : call === "PARTIAL" ? "Partial" : "Absent"}
    </span>
  );
}

export function Modal({
  open,
  onClose,
  title,
  children,
  wide,
}: {
  open: boolean;
  onClose: () => void;
  title: ReactNode;
  children: ReactNode;
  wide?: boolean;
}) {
  useEffect(() => {
    if (!open) return;
    const onKey = (e: KeyboardEvent) => {
      if (e.key === "Escape") onClose();
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [open, onClose]);
  if (!open) return null;
  return (
    <div
      className="fixed inset-0 z-50 flex items-start justify-center bg-zinc-900/40 p-4 sm:p-8 overflow-y-auto"
      onMouseDown={(e) => {
        if (e.target === e.currentTarget) onClose();
      }}
    >
      <div
        className={`bg-white rounded-xl shadow-xl border border-zinc-200 w-full dark:bg-zinc-900 dark:border-zinc-800 ${wide ? "max-w-5xl" : "max-w-2xl"} my-auto`}
      >
        <div className="flex items-center justify-between px-6 py-4 border-b border-zinc-200 dark:border-zinc-800">
          <h2 className="text-lg font-semibold">{title}</h2>
          <button
            onClick={onClose}
            className="w-11 h-11 grid place-items-center rounded-md text-zinc-500 hover:text-zinc-900 hover:bg-zinc-100 dark:text-zinc-400 dark:hover:text-zinc-100 dark:hover:bg-zinc-800"
            aria-label="Close"
          >
            <svg width="16" height="16" viewBox="0 0 16 16" fill="none">
              <path
                d="M3 3l10 10M13 3L3 13"
                stroke="currentColor"
                strokeWidth="1.8"
                strokeLinecap="round"
              />
            </svg>
          </button>
        </div>
        <div className="p-6">{children}</div>
      </div>
    </div>
  );
}

export function Spinner({ className = "" }: { className?: string }) {
  return (
    <svg
      className={`animate-spin ${className}`}
      width="20"
      height="20"
      viewBox="0 0 24 24"
      fill="none"
      aria-label="Loading"
    >
      <circle
        className="opacity-20"
        cx="12"
        cy="12"
        r="10"
        stroke="currentColor"
        strokeWidth="4"
      />
      <path
        className="opacity-90"
        fill="currentColor"
        d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z"
      />
    </svg>
  );
}

export function ErrorBox({ message }: { message: string }) {
  return (
    <div className="rounded-lg border border-red-200 bg-red-50 text-red-800 px-4 py-3 text-[15px] dark:border-red-900 dark:bg-red-950/40 dark:text-red-300">
      <span className="font-semibold">Something went wrong. </span>
      {message}
    </div>
  );
}

export function InfoIcon({ text }: { text: string }) {
  return (
    <span className="relative group inline-flex">
      <svg
        width="16"
        height="16"
        viewBox="0 0 16 16"
        className="text-zinc-400 dark:text-zinc-500"
        aria-label={text}
      >
        <circle cx="8" cy="8" r="7" fill="none" stroke="currentColor" strokeWidth="1.4" />
        <rect x="7.3" y="6.5" width="1.4" height="4.6" rx="0.7" fill="currentColor" />
        <circle cx="8" cy="4.4" r="0.9" fill="currentColor" />
      </svg>
      <span className="pointer-events-none absolute left-1/2 -translate-x-1/2 bottom-full mb-2 w-64 rounded-md bg-zinc-900 text-white text-xs leading-relaxed px-3 py-2 opacity-0 group-hover:opacity-100 transition-opacity z-20 dark:bg-zinc-100 dark:text-zinc-900">
        {text}
      </span>
    </span>
  );
}

/** Drag and drop upload area with click to browse. */
export function DropZone({
  onFiles,
  multiple = false,
  accept = ".fasta,.fa,.fna,.fsa",
  hint,
  compact,
}: {
  onFiles: (files: File[]) => void;
  multiple?: boolean;
  accept?: string;
  hint: string;
  compact?: boolean;
}) {
  const inputRef = useRef<HTMLInputElement>(null);
  const [dragging, setDragging] = useState(false);
  function handleDrop(e: React.DragEvent) {
    e.preventDefault();
    setDragging(false);
    const files = Array.from(e.dataTransfer.files);
    if (files.length) onFiles(files);
  }
  return (
    <div
      onDragOver={(e) => {
        e.preventDefault();
        setDragging(true);
      }}
      onDragLeave={() => setDragging(false)}
      onDrop={handleDrop}
      onClick={() => inputRef.current?.click()}
      role="button"
      tabIndex={0}
      onKeyDown={(e) => {
        if (e.key === "Enter" || e.key === " ") inputRef.current?.click();
      }}
      className={`cursor-pointer rounded-xl border-2 border-dashed transition-colors text-center ${
        compact ? "p-4" : "p-8"
      } ${
        dragging
          ? "border-zinc-900 bg-zinc-100 dark:border-zinc-100 dark:bg-zinc-800"
          : "border-zinc-300 bg-white hover:border-zinc-400 hover:bg-zinc-50 dark:border-zinc-700 dark:bg-zinc-900 dark:hover:border-zinc-600 dark:hover:bg-zinc-800/60"
      }`}
    >
      <input
        ref={inputRef}
        type="file"
        accept={accept}
        multiple={multiple}
        className="hidden"
        onChange={(e) => {
          const files = Array.from(e.target.files ?? []);
          if (files.length) onFiles(files);
          e.target.value = "";
        }}
      />
      <svg
        width="28"
        height="28"
        viewBox="0 0 24 24"
        fill="none"
        className={`mx-auto ${compact ? "mb-2" : "mb-3"} text-zinc-400 dark:text-zinc-500`}
      >
        <path
          d="M12 16V4m0 0l-4 4m4-4l4 4M4 20h16"
          stroke="currentColor"
          strokeWidth="1.8"
          strokeLinecap="round"
          strokeLinejoin="round"
        />
      </svg>
      <p className="text-[15px] text-zinc-700 dark:text-zinc-300">
        Drag files here, or <span className="text-zinc-900 underline dark:text-zinc-100">browse</span>
      </p>
      <p className="text-sm text-zinc-500 mt-1 dark:text-zinc-500">{hint}</p>
    </div>
  );
}

import { useState } from "react";

export function Tabs({
  tabs,
  active,
  onChange,
}: {
  tabs: { key: string; label: string; disabled?: boolean }[];
  active: string;
  onChange: (key: string) => void;
}) {
  return (
    <div className="flex gap-1 border-b border-zinc-200 dark:border-zinc-800" role="tablist">
      {tabs.map((t) => (
        <button
          key={t.key}
          role="tab"
          aria-selected={active === t.key}
          disabled={t.disabled}
          onClick={() => onChange(t.key)}
          className={`px-5 h-12 text-[15px] font-medium rounded-t-lg border-b-2 -mb-px transition-colors disabled:text-zinc-300 disabled:cursor-not-allowed dark:disabled:text-zinc-700 ${
            active === t.key
              ? "border-zinc-900 text-zinc-900 dark:border-zinc-100 dark:text-zinc-100"
              : "border-transparent text-zinc-500 hover:text-zinc-900 dark:text-zinc-400 dark:hover:text-zinc-100"
          }`}
        >
          {t.label}
        </button>
      ))}
    </div>
  );
}

export function EmptyState({
  title,
  hint,
  action,
}: {
  title: string;
  hint?: string;
  action?: ReactNode;
}) {
  return (
    <div className="text-center py-16 px-6">
      <p className="text-lg font-medium text-zinc-700 dark:text-zinc-300">{title}</p>
      {hint && <p className="text-[15px] text-zinc-500 mt-2 max-w-md mx-auto dark:text-zinc-400">{hint}</p>}
      {action && <div className="mt-6 flex justify-center">{action}</div>}
    </div>
  );
}
