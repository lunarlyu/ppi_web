import { cp, mkdir, rm, writeFile } from "node:fs/promises";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const root = join(dirname(fileURLToPath(import.meta.url)), "..");
const dist = join(root, "dist");

await rm(dist, { force: true, recursive: true });
await mkdir(join(dist, "client"), { recursive: true });
await mkdir(join(dist, "server"), { recursive: true });
await mkdir(join(dist, ".openai"), { recursive: true });

await cp(join(root, "index.html"), join(dist, "client", "index.html"));
await cp(join(root, "assets"), join(dist, "client", "assets"), {
  recursive: true,
});
await cp(join(root, "data"), join(dist, "client", "data"), {
  recursive: true,
});
await cp(join(root, "public"), join(dist, "client"), {
  recursive: true,
});
await cp(join(root, ".openai", "hosting.json"), join(dist, ".openai", "hosting.json"));

await writeFile(
  join(dist, "server", "index.js"),
  `export default {
  async fetch(request, env) {
    const url = new URL(request.url);
    const response = await env.ASSETS.fetch(request);

    if (response.status !== 404 || url.pathname.includes(".")) {
      return response;
    }

    return env.ASSETS.fetch(new URL("/", request.url));
  },
};
`,
);
