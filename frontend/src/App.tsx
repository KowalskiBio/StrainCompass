import { Link, Navigate, Route, Routes, useLocation } from "react-router-dom";
import ProjectsPage from "./pages/ProjectsPage";
import ProjectPage from "./pages/ProjectPage";
import SettingsPage from "./pages/SettingsPage";

export default function App() {
  const location = useLocation();
  const onSettings = location.pathname.startsWith("/settings");
  return (
    <div className="min-h-screen flex flex-col">
      <header className="bg-white border-b border-zinc-200">
        <div className="max-w-7xl mx-auto px-6 h-16 flex items-center justify-between">
          <Link to="/" className="flex items-center gap-3 group">
            <span className="w-8 h-8 rounded-md bg-zinc-900 text-white grid place-items-center font-bold text-sm">
              b
            </span>
            <span className="text-lg font-semibold tracking-tight">bactiment</span>
          </Link>
          <nav className="flex items-center gap-1">
            <Link
              to="/"
              className={`px-4 h-11 inline-flex items-center rounded-md text-[15px] font-medium transition-colors ${
                !onSettings
                  ? "text-zinc-900"
                  : "text-zinc-500 hover:text-zinc-900 hover:bg-zinc-100"
              }`}
            >
              Projects
            </Link>
            <Link
              to="/settings"
              className={`px-4 h-11 inline-flex items-center rounded-md text-[15px] font-medium transition-colors ${
                onSettings
                  ? "text-zinc-900"
                  : "text-zinc-500 hover:text-zinc-900 hover:bg-zinc-100"
              }`}
            >
              Settings
            </Link>
          </nav>
        </div>
      </header>
      <main className="flex-1">
        <Routes>
          <Route path="/" element={<ProjectsPage />} />
          <Route path="/projects/:id" element={<ProjectPage />} />
          <Route path="/settings" element={<SettingsPage />} />
          <Route path="*" element={<Navigate to="/" replace />} />
        </Routes>
      </main>
      <footer className="border-t border-zinc-200 bg-white">
        <div className="max-w-7xl mx-auto px-6 py-3 text-sm text-zinc-400">
          bactiment: bacterial genome comparison workbench
        </div>
      </footer>
    </div>
  );
}
