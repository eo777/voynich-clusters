// Node's ESM loader requires file URLs for Windows paths. WebR 0.6.0 passes
// an absolute Windows path when importing its R runtime in a worker.
import { pathToFileURL } from 'node:url';
export async function resolve(specifier, context, nextResolve) {
  const resolved = /^[A-Za-z]:[\\/]/.test(specifier)
    ? pathToFileURL(specifier).href : specifier;
  return nextResolve(resolved, context);
}
