# RR5e.0 resource lifecycle audit

Date: 2026-08-11

The stale failure was not reproducible after removal of the ACL-inaccessible pytest roots and execution in the normal Windows test context.

Validation evidence:

- first ordinary full suite: `182 passed in 42.04s`;
- after the first suite: zero Python/pytest processes;
- the full basetemp and cache were recursively deletable;
- both directories were recreated and deleted again successfully;
- second ordinary full suite from that clean state: `182 passed in 42.13s`;
- after the second suite: zero Python/pytest processes;
- an actual generated `manifest.json` was renamed to a probe name and back successfully;
- the complete basetemp and cache were then recursively deleted successfully.

These checks establish that process pools shut down, workers joined, generated manifests/checkpoints were closed, and the test task bundles no longer held rename/delete-blocking handles. The existing resumability tests also passed in both suites, including interruption/resume equivalence, fingerprint mismatch refusal, corrupt/incomplete bundle quarantine, and terminal-set accounting.

No reproducible RR5d resource-lifecycle bug was established. Consequently no scientific source, engineering source, or test code was changed, and no regression test was added merely to encode an external ACL condition.
