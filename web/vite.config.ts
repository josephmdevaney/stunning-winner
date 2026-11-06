import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// The API target is overridable so a second checkout (a git worktree) can run its own backend on another port.
const apiTarget = process.env.VITE_API_TARGET ?? 'http://127.0.0.1:8000'

// https://vite.dev/config/
export default defineConfig({
  plugins: [react()],
  server: {
    proxy: {
      '/api': apiTarget,
    },
  },
})
