# obflow bot

ChatOps for maintaining the benchmark plan (`benchmark_conda.yaml`). Comment a
command on **any PR or issue** and the bot rewrites the plan's module commits for
you — no hand-editing 20 YAML fields.

Implemented by [`obflow-bot.yml`](./obflow-bot.yml); the actual rewrite is done
by [`obflow`](https://github.com/btraven00/obflow) (`plan pin` / `plan track`).

## Commands

| Comment | What it does |
|---|---|
| `/pin-commits` | **Freeze.** Dereference each module's tracked branch to a concrete commit SHA. |
| `/pin-commits <ref>` | Freeze every module to `<ref>`'s SHA instead of its own branch. |
| `/track-branch <branch>` | **Thaw.** Point each module's commit at `<branch>` — for every module whose repo actually has that branch. |

`pin` and `track` are inverses: `pin` moves branch → SHA, `track` moves SHA (or
another branch) → branch.

### Where the change lands

- **On a PR** → committed straight onto that PR's branch (rides along with the
  work under review).
- **On an issue** → opens/updates a dedicated PR (branch `obflow/plan-update`).

The bot 👀-reacts when it picks up the command, then 👍/👎-reacts and replies
with a per-module report when done.

## Behavior worth knowing

- **A frozen module is never silently re-pinned.** `/pin-commits` skips any
  module already at a SHA (reported `already-pinned`) — it only resolves modules
  currently tracking a branch. Re-running it is a no-op.
- **`/track-branch` only touches modules that have the branch.** Repos without it
  are left exactly as they are (reported `not-found`).
- **SHAs are abbreviated** to 7 chars, matching the plan's existing convention.
- **Minimal diffs.** Only `commit:` lines change; indentation, comments, and
  blank lines are preserved byte-for-byte.

Each module gets a **status**: `pinned`, `already-pinned`, `tracking`,
`unchanged`, `not-found`, or `error` — shown in the reply/PR table (the bot reads
obflow's `--json` output).

## Who can run it

Only commenters with write access (`OWNER` / `MEMBER` / `COLLABORATOR`). Others
are ignored.

## One-time setup (for a fork / new repo)

1. **Default branch.** This workflow must live on the repo's **default branch** —
   `issue_comment` events only trigger workflows from there.
2. **Actions settings** (Settings → Actions → General → bottom):
   - Workflow permissions → **Read and write permissions**
   - ☑ **Allow GitHub Actions to create and approve pull requests**
3. **Enable Actions and Issues** on the fork if they aren't already.
4. **obflow binary.** The bot downloads `obflow` from its `nightly` release
   (`OBFLOW_URL` in the workflow). The `--remote`/`track`/`--abbrev` features must
   be present in that build — rebuild the release after updating obflow.
5. **Private modules** (optional). If any module repo is private, add a
   `MODULES_TOKEN` secret (a PAT with read access); it's wired into git so
   `ls-remote` can reach them. Public modules need nothing.
