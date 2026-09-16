import { Link, Navigate, Route, Routes, useLocation } from "react-router-dom";
import ProjectsPage from "./pages/ProjectsPage";
import ProjectPage from "./pages/ProjectPage";
import SettingsPage from "./pages/SettingsPage";
import { useTheme } from "./theme";

export default function App() {
  const location = useLocation();
  const onSettings = location.pathname.startsWith("/settings");
  return (
    <div className="min-h-screen flex flex-col">
      <header className="bg-white border-b border-zinc-200 dark:bg-zinc-900 dark:border-zinc-800">
        <div className="max-w-[1600px] mx-auto px-6 h-16 flex items-center justify-between">
          <Link to="/" className="flex items-center gap-3 group">
            <span className="w-8 h-8 rounded-md bg-zinc-900 text-white grid place-items-center font-bold text-sm dark:bg-zinc-100 dark:text-zinc-900">
              b
            </span>
            <span className="text-lg font-semibold tracking-tight">straincompass</span>
          </Link>
          <div className="flex items-center gap-1">
            <nav className="flex items-center gap-1">
              <Link
                to="/"
                className={`px-4 h-11 inline-flex items-center rounded-md text-[15px] font-medium transition-colors ${
                  !onSettings
                    ? "text-zinc-900 dark:text-zinc-100"
                    : "text-zinc-500 hover:text-zinc-900 hover:bg-zinc-100 dark:text-zinc-400 dark:hover:text-zinc-100 dark:hover:bg-zinc-800"
                }`}
              >
                Projects
              </Link>
              <Link
                to="/settings"
                className={`px-4 h-11 inline-flex items-center rounded-md text-[15px] font-medium transition-colors ${
                  onSettings
                    ? "text-zinc-900 dark:text-zinc-100"
                    : "text-zinc-500 hover:text-zinc-900 hover:bg-zinc-100 dark:text-zinc-400 dark:hover:text-zinc-100 dark:hover:bg-zinc-800"
                }`}
              >
                Settings
              </Link>
            </nav>
            <ThemeToggle />
          </div>
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
      <footer className="border-t border-zinc-200 bg-white dark:border-zinc-800 dark:bg-zinc-900">
        <div className="max-w-[1600px] mx-auto px-6 py-3 text-sm text-zinc-400 dark:text-zinc-500">
          straincompass: bacterial genome comparison workbench
        </div>
      </footer>
    </div>
  );
}

function ThemeToggle() {
  const { theme, toggle } = useTheme();
  return (
    <button
      onClick={toggle}
      aria-label={theme === "dark" ? "Switch to light mode" : "Switch to dark mode"}
      title={theme === "dark" ? "Switch to light mode" : "Switch to dark mode"}
      className="w-11 h-11 grid place-items-center rounded-md text-zinc-500 hover:text-zinc-900 hover:bg-zinc-100 dark:text-zinc-400 dark:hover:text-zinc-100 dark:hover:bg-zinc-800"
    >
      {theme === "dark" ? (
        <svg width="18" height="18" viewBox="0 0 16 16" fill="none" aria-hidden>
          <path
            d="M13.5 9.3A5.6 5.6 0 016.7 2.5a5.6 5.6 0 106.8 6.8z"
            stroke="currentColor"
            strokeWidth="1.4"
            strokeLinejoin="round"
          />
        </svg>
      ) : (
        <svg width="18" height="18" viewBox="0 0 16 16" fill="none" aria-hidden>
          <circle cx="8" cy="8" r="3.2" stroke="currentColor" strokeWidth="1.4" />
          <path
            d="M8 1v1.5M8 13.5V15M15 8h-1.5M2.5 8H1M12.7 3.3l-1.1 1.1M4.4 11.6l-1.1 1.1M12.7 12.7l-1.1-1.1M4.4 4.4L3.3 3.3"
            stroke="currentColor"
            strokeWidth="1.4"
            strokeLinecap="round"
          />
        </svg>
      )}
    </button>
  );
}
