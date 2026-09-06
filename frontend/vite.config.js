import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import path from "path";

// https://vite.dev/config/
export default defineConfig({
  plugins: [react()],
  resolve: {
    alias: {
      // ensure imports like "plotly.js/dist/plotly" resolve to the actual minified file
      "plotly.js/dist/plotly": path.resolve(
        __dirname,
        "node_modules/plotly.js-dist-min/plotly.min.js"
      ),
      "plotly.js": path.resolve(
        __dirname,
        "node_modules/plotly.js-dist-min/plotly.min.js"
      ),
    },
  },
  server: {
    proxy: {
      "/api": {
        target: "http://localhost:8001",
        changeOrigin: true,
      },
    },
  },
});
