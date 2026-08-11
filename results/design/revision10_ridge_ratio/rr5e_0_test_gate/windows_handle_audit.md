# RR5e.0 Windows handle and permission audit

Date: 2026-08-11

## Process audit

Before cleanup, ordinary `Get-Process` and elevated `Get-CimInstance Win32_Process` checks found no Python, pytest, multiprocessing worker, or process-pool process. No process was terminated. The same elevated process audit after each successful full suite found zero Python/pytest processes.

## Filesystem audit

Both inaccessible paths existed and were directories:

- `results/pytest-parallel-audit`, created 2026-08-11 08:46:40 and last written 2026-08-11 08:46:54;
- `.pytest_cache`, created and last written initially on 2026-08-09.

Neither root had a read-only attribute. Elevated enumeration showed pytest node directories, pytest `current` reparse points, smoke-test outputs, and cache-provider files only. Git tracks neither path, and `.gitignore` excludes `results/pytest-*` and `.pytest_cache*`.

`icacls` reported the following effective ACL shape on both stale roots:

```text
NT AUTHORITY\SYSTEM:(OI)(CI)(F)
BUILTIN\Administrators:(OI)(CI)(F)
OWNER RIGHTS:(OI)(CI)(F)
```

The ordinary repository process could not list either directory or even read its ACL. Elevated access could enumerate and delete both immediately. No sharing-violation/open-handle error occurred.

## Classification

Root cause: **PERMISSION/ACL plus STALE TEST ARTIFACT**.

Not found: open-handle lock, read-only attribute, live orphan worker, or scientific-file involvement. The evidence does not support an RR5d process-pool or file-handle leak.
