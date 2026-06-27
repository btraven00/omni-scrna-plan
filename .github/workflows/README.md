# obflow bot

ChatOps for maintaining the benchmark plan and its modules. Comment a command on
a PR or issue and the bot does the multi-repo edit for you — no hand-editing
20 YAML fields across repos.

Implemented by [`obflow-bot.yml`](./obflow-bot.yml) (pin/track) and
[`obflow-retarget.yml`](./obflow-retarget.yml) (retarget); the actual edits are
done by [`obflow`](https://github.com/btraven00/obflow).

## Commands

| Comment | What it does |
|---|---|
| `/pin-commits` | **Freeze.** Dereference each module's tracked branch to a concrete commit SHA, in the plan. |
| `/pin-commits <ref>` | Freeze every module to `<ref>`'s SHA instead of its own branch. |
| `/track-branch <branch>` | **Thaw.** Point each module's commit at `<branch>` (in the plan) — for every module whose repo has that branch. |
| `/retarget-modules [target] [limit]` | **Land.** In each module repo, open a PR setting *this plan's* benchmark `ref` (in its own `omnibenchmark.yaml`) to `target` (default `main`). `limit` caps how many PRs to open (debug). |

`pin` and `track` are inverses on the **plan** (`benchmark_conda.yaml`): `pin`
moves branch → SHA, `track` moves SHA (or another branch) → branch.
`retarget-modules` is the post-merge cleanup on the **modules' own**
`omnibenchmark.yaml` files (a different file, in a different repo).

### Where the change lands

- `/pin-commits`, `/track-branch`:
  - **On a PR** → committed straight onto that PR's branch.
  - **On an issue** → opens/updates a dedicated PR (branch `obflow/plan-update`).
- `/retarget-modules` → opens one PR **per module repo** (branch
  `obflow/retarget-<target>`), and replies with a repo-by-repo table.

The bot 👀-reacts when it picks up the command, then 👍/👎-reacts and replies
with a per-module report when done.

## Behavior worth knowing

- **A frozen module is never silently re-pinned.** `/pin-commits` skips any
  module already at a SHA (reported `already-pinned`) — it only resolves modules
  currently tracking a branch. Re-running it is a no-op.
- **`/track-branch` only touches modules that have the branch.** Repos without it
  are left exactly as they are (reported `not-found`).
- **SHAs are abbreviated** to 7 chars, matching the plan's existing convention.
- **`/retarget-modules` only touches this plan's benchmark entry** in each module's
  `omnibenchmark.yaml` — `templates.ref` and any other benchmarks are left alone.
  Modules with no matching entry are reported `no-match` (no PR). Idempotent.
- **Minimal diffs.** Only the relevant line(s) change; indentation, comments, and
  blank lines are preserved byte-for-byte.

Each module gets a **status** — `pinned`, `already-pinned`, `tracking`,
`unchanged`, `not-found`, `retargeted`, `ref-differs`, `no-match`, or `error` —
shown in the reply/PR table (the bot reads obflow's `--json` output).

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
5. **Private modules** (optional, pin/track). If any module repo is private, add a
   `MODULES_TOKEN` secret (a PAT with read access); it's wired into git so
   `ls-remote` can reach them. Public modules need nothing.
6. **The omni-clanker GitHub App (required for `/retarget-modules`).** That
   command writes to the *module* repos, which the default `GITHUB_TOKEN` cannot
   reach. See [Creating the omni-clanker App](#creating-the-omni-clanker-app)
   below; the workflow mints a short-lived token from it at runtime. Without the
   App secrets, `/retarget-modules` replies with an error and does nothing.

## Creating the omni-clanker App

`/retarget-modules` opens PRs in the *module* repos. The plan repo's
`GITHUB_TOKEN` can't do that, so we use a GitHub App owned by the org as a scoped,
expiring credential. The App is just an identity + key — **no server, no webhook,
no hosted code**; the bot logic stays in this workflow.

**1. Register it** (you must be an `omni-scrna` org owner):
GitHub → `omni-scrna` org → Settings → Developer settings → GitHub Apps →
**New GitHub App**.
- **Name:** `omni-clanker`
- **Homepage URL:** anything (e.g. this repo's URL)
- **Webhook:** **uncheck "Active"** (not needed)
- **Repository permissions:** *Contents* → Read and write, *Pull requests* →
  Read and write. Leave everything else "No access".
- **Where can this app be installed:** "Only on this account" (keep it private).
- Create the App, then note the **App ID**, and **Generate a private key**
  (downloads a `.pem`).

**2. Install it:** on the App's page → **Install App** → `omni-scrna` → select the
module repos (or "All repositories").

**3. Add the secrets** to this repo (Settings → Secrets and variables → Actions),
and to your fork if testing there:
- `OMNI_CLANKER_APP_ID` — the numeric App ID
- `OMNI_CLANKER_PRIVATE_KEY` — the entire contents of the `.pem` file

The workflow's `actions/create-github-app-token` step exchanges these for a
~1-hour installation token scoped to the App's permissions and installed repos.
PRs it opens are authored by `omni-clanker[bot]`.
