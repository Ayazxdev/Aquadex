function resolveApiBase() {
  const envBase = import.meta.env.VITE_API_BASE;
  if (!envBase) {
    return "http://localhost:8001/api";
  }

  // Remove trailing slashes and spaces
  let base = envBase.trim().replace(/\/+$/, "");

  // If user passed only domain (e.g. https://xxx.onrender.com), append /api
  if (!base.endsWith("/api")) {
    base = `${base}/api`;
  }

  return base;
}

export const API_BASE = resolveApiBase();
console.log("[Aquadex] Configured API_BASE:", API_BASE);
