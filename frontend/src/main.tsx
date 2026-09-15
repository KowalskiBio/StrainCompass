import React from "react";
import { createRoot } from "react-dom/client";
import { BrowserRouter, Route, Routes } from "react-router-dom";
import App from "./App";
import "./index.css";

/** A crash shows this card instead of a white screen. */
class ErrorBoundary extends React.Component<
  { children: React.ReactNode },
  { error: Error | null }
> {
  state = { error: null as Error | null };

  static getDerivedStateFromError(error: Error) {
    return { error };
  }

  render() {
    if (this.state.error) {
      return (
        <div className="min-h-screen flex items-center justify-center p-6">
          <div className="max-w-lg w-full border border-red-200 bg-red-50 rounded-xl p-6">
            <h1 className="text-lg font-semibold text-red-900">
              Something went wrong
            </h1>
            <p className="mt-2 text-[15px] text-red-800">
              The page hit an unexpected error. Reloading usually fixes it;
              your data on the server is safe.
            </p>
            <pre className="mt-3 text-xs text-red-700 bg-white/60 rounded-lg p-3 overflow-auto max-h-40 font-mono">
              {this.state.error.message}
            </pre>
            <div className="mt-4 flex gap-3">
              <button
                className="h-11 px-4 rounded-lg bg-zinc-900 text-white text-[15px] hover:bg-zinc-700"
                onClick={() => window.location.reload()}
              >
                Reload the page
              </button>
              <button
                className="h-11 px-4 rounded-lg border border-red-300 bg-white text-red-700 text-[15px] hover:bg-red-100"
                onClick={() => this.setState({ error: null })}
              >
                Try again
              </button>
            </div>
          </div>
        </div>
      );
    }
    return this.props.children;
  }
}

createRoot(document.getElementById("root")!).render(
  <React.StrictMode>
    <ErrorBoundary>
      <BrowserRouter>
        <Routes>
          {/* App renders the header/footer and resolves the real routes
              (/, /projects/:id, /settings) in its own Routes */}
          <Route path="/*" element={<App />} />
        </Routes>
      </BrowserRouter>
    </ErrorBoundary>
  </React.StrictMode>,
);
