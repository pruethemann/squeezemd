# CLAUDE.md — Project Guidelines for Claude Code

This file is read automatically by Claude Code at the start of every session.
Follow all instructions in this file unless the user explicitly overrides them in the current session.

---

## 1. Mode Overview

Claude Code operates in two distinct modes. Always be explicit about which mode is active.

| Mode | Activated by | File writes allowed? |
|---|---|---|
| **Plan Mode** | `Shift+Tab` twice, or `/plan` | ❌ Read-only (except `claudecode-instructions/`) |
| **Execution Mode** | Default / after plan approval | ✅ Full access |

---

## 2. Plan Mode Guidelines

> **Plan Mode is read-only.** Do not create, edit, or delete any file outside of `claudecode-instructions/`.

### 2a. Model Selection — No Restriction, But Be Smart

- You may use any available model in Plan Mode (Opus, Sonnet, or Haiku).
- **Prefer not to call Opus unless the task genuinely requires it** (see model guide below).
- If you are unsure, default to Sonnet.

**Model Selection Guide:**

| Task Complexity | Recommended Model | Examples |
|---|---|---|
| Simple / well-scoped | **Haiku** | Single-file lookup, explaining a function, quick grep |
| Moderate / multi-file | **Sonnet** *(default)* | Planning a feature, understanding module interactions |
| High complexity / architectural | **Opus** | Refactoring across many modules, resolving deep ambiguity, novel architecture decisions |

Always state your model recommendation at the end of a Plan Mode response:
> 💡 *Recommended model for execution: **Sonnet** — this task involves N files with clear scope.*

### 2b. Knowledge Graph Lookup (Graphify)

- If the task requires looking up **more than 3 different files**, consult the knowledge graph in `graphify-out/` **before** reading individual files.
- The knowledge graph provides a structural overview of the codebase and reduces context usage.
- Usage pattern:
  1. Check `graphify-out/` for the relevant module or entity.
  2. Use the graph to identify the minimal set of files to read.
  3. Only then open specific files for detail.

> ⚠️ If `graphify-out/` is empty or stale, note this in your plan and suggest the user run `graphify --update` before execution.

---

## 3. Execution Mode Guidelines

> **Execution Mode makes real changes.** Proceed carefully and verify after every run.

### 3a. Opus Warning

If Opus is selected or requested for execution:

1. **Display a warning before proceeding:**
   ```
   ⚠️  WARNING: You are about to run in Execution Mode with Opus.
   Opus is the most capable and most expensive model.
   Do you want to continue? (yes / no)
   ```
2. **Wait for explicit confirmation** (`yes`) before taking any action.
3. If the user responds with anything other than `yes`, abort and suggest Sonnet as an alternative.

### 3b. Post-Execution Git Diff Check

After **every execution run**, immediately run:

```bash
git diff HEAD --shortstat
```

Then evaluate the output:

| Condition | Action |
|---|---|
| ≤ 5 files changed **and** ≤ 200 lines changed | ✅ Continue normally |
| > 5 files created/modified **or** > 200 lines changed | 🔶 Trigger commit recommendation (see below) |

**Commit recommendation message (display when threshold is hit):**

```
🔶 CHECKPOINT RECOMMENDED

Since the last commit:
  • Files changed: X  (threshold: 5)
  • Lines changed: Y  (threshold: 200)

Recommended actions:
  1. Review the diff:    git diff HEAD
  2. Commit the work:    git add -A && git commit -m "your message"
  3. Update the graph:   graphify --update

Would you like help writing the commit message?
```

> **Rationale:** Keeping commits small makes it easier to revert if something goes wrong, and keeps the knowledge graph in sync with the actual codebase.

---

## 4. Code Organisation Guidelines

**Store large generated data in `data/`.** Any large data files directly produced by scripts must be written to `data/`, not left alongside source files or in the project root.

**Keep scripts short and modular.** Avoid generating scripts longer than ~2,000 lines. Prefer implementing logic as reusable functions or modules in `src/`, then calling them from `scripts/` or `analysis/` files. This keeps individual files reviewable, testable, and easy to diff.

---

## 5. Session Management Guidelines

**One task per session.** Close a session and open a new one when a task is done. Resist the temptation to continue into unrelated work in the same session. Your knowledge graph (`graphify-out/`) and Claude memory (`.docs/`) are designed to carry context forward — you should not need the chat history.

**Compact when a session grows long.** If a single session has gone through many exchanges, use `/compact` to summarise the conversation before continuing. Signs you should compact: responses feel slower, Claude starts repeating itself, or the task has gone through more than one major pivot.

**Finish one thing before starting another.** A session scoped to a single task is easier to review, easier to roll back, and produces a cleaner instruction file. If you find yourself doing two unrelated things in one session, finish the first, commit, and open a fresh session for the second.

**Write down what matters.** After a session ends, any important decisions or context should be captured in `.docs/` or a `claudecode-instructions/` file — not left in chat history. If Claude needs to know something next time, write it down explicitly. Chat history is ephemeral; the file system is not.

**Commit and update the graph before closing.** Before ending a session, check the git diff and commit anything meaningful. Then run `graphify --update` if files were added or restructured. This ensures the next session starts from a clean, accurate baseline.

---

## 6. Folder Conventions (Quick Reference)

| Folder | Who writes to it | Notes |
|---|---|---|
| `data/` | Claude (execution) + scripts | Large generated data files |
| `src/` | Claude (execution) | Reusable modules and functions |
| `scripts/` | Claude (execution) | Entry-point scripts; call into `src/` |
| `analysis/` | Claude (execution) | Analysis files; call into `src/` |
| `skills/` | Humans + Claude (execution) | Reusable Claude Code skill `.md` files |
| `claudecode-instructions/` | Claude (plan mode) + Humans | Task instructions and plans |
| `.docs/` | Claude (execution) | Persistent memory; do not edit manually |
| `graphify-out/` | Graphify tool | Read-only for Claude; update via `graphify --update` |

---

## 7. General Principles

- **Plan before you act.** Use Plan Mode for any task spanning more than one file or module.
- **Small, reviewable changes.** Prefer many small commits over one large one.
- **Flag uncertainty early.** If requirements are ambiguous, say so in the plan — do not guess during execution.
- **Prefer reversibility.** When in doubt, choose the approach that is easiest to undo.
- **Keep this file updated.** If the team agrees on a new convention, add it here.