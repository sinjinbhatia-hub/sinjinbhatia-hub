# 30-Day Agentic Build Plan — Wardrobe Try-On PWA

Companion to `feasibility-report.md` (Deliverable 4 in executable form). This file is written
to be handed to Claude Code agents task-by-task. Each task states: goal, definition of done,
and how the founder verifies it **from a phone** (no laptop at any step).

---

## 0. Stack decision (locked)

**Next.js 15 (App Router, TypeScript) on Vercel Hobby + Supabase Free (Postgres, Auth, Storage) + GitHub Actions.**
No FastAPI, no Railway in v1.

Why this and not the stack you already know end-to-end:

- **One deployable.** Every PR gets a Vercel preview URL you open in mobile Safari. A separate
  FastAPI service doubles the deploy surfaces, secrets, CORS config, and failure modes — and
  there is no genuinely free, cold-start-free Python host left (Railway = one-time $5 trial;
  Render free = ~50s cold starts, which kills a demo).
- **Server code still exists** — Next.js API routes (Node) hold the API keys and call
  Gemini / fal. You never expose keys to the browser.
- **You know React**; Next is React plus file-routing. The agents writing the code are highly
  fluent in Next.js — on a phone, you want the stack the *agents* are best at, because you
  can't drop into a debugger to save them.
- **Python escape hatch:** if a self-hosted try-on model is ever needed, it runs on Modal
  (deployed from GitHub Actions, verified via hosted endpoint) — not on a server you babysit.

### Architecture (generation is async, always)

```
[PWA / mobile Safari]
   │  upload photo (Supabase Storage, signed URL)
   ▼
[Next.js API route]  ──create──▶  jobs table (Supabase)
   │                                   ▲
   │  call provider (Gemini image API / fal queue)   │
   ▼                                   │
[provider]  ──webhook or server poll──┘
   │
   ▼
result image → Supabase Storage → job row updated → client polls /api/jobs/:id
```

Rules that keep this debuggable from a phone:

1. **Every generation is a `jobs` row** (id, user, type, provider, status, input refs, output
   ref, error, timing, cost estimate). The Supabase table editor is your ops dashboard.
2. **Client polls; server never blocks >10s.** Works regardless of Vercel function limits.
3. **Provider abstraction from day one**: `lib/providers/{gemini,fal,fashn}.ts` behind one
   `generateTryOn(input): jobId` interface, provider chosen per-job by a config flag — so A/B
   between engines is a env-var flip, not a refactor.
4. **Canonical avatar invariant:** all try-ons are a SINGLE edit from the stored canonical
   avatar image. Never edit an edited image (identity drift compounds). Multi-garment outfits
   = one multi-image call, not sequential dressing.

## 1. Repo structure

```
/                      Next.js app (app router)
  app/
    (marketing)/       landing page
    onboard/           photo upload → avatar creation → consent
    closet/            garment grid, add-garment flow
    dressing-room/     pick garment(s) → generate → result viewer
    shop/              recommendation rail + affiliate links
    api/
      jobs/            create/status
      webhooks/fal/    provider callbacks
      health/          selftest endpoint (env, db, provider key checks)
  lib/
    providers/         gemini.ts, fal.ts, fashn.ts (one interface)
    style/             art-direction prompt templates + presets (VERSIONED)
    catalog/           attribute extraction, bg removal client hook
  supabase/
    migrations/        SQL migrations (single-writer — see parallelization)
    seed/              demo user + demo closet script
  e2e/                 Playwright smoke tests (run in CI against preview URL)
  .github/workflows/
    ci.yml             typecheck, lint, build, e2e vs preview, screenshot artifacts
    seed.yml           workflow_dispatch: seed demo data
    lookbook.yml       workflow_dispatch: batch pre-generate outfit images
    keepalive.yml      weekly ping so Supabase free project never pauses
  docs/                this plan, the feasibility report, decision log
```

## 2. Phone-verification toolkit (set up in week 1, use forever)

| What | How you check it from the phone |
|---|---|
| Every code change | Vercel preview URL on the PR (mobile Safari) |
| Visual regressions | Playwright screenshots in CI, committed to the PR branch under `e2e/screens/` — GitHub mobile renders them in the diff |
| Backend state | Supabase dashboard (tables, storage, auth, logs) in mobile browser |
| Env/config sanity | `/api/health` returns JSON: env vars present? DB reachable? provider key valid (1¢ test call)? |
| Batch operations | GitHub Actions `workflow_dispatch` buttons — pressable from the GitHub mobile app (seed data, regenerate lookbook, run smoke suite) |
| Generation failures | `jobs` table: status, error text, provider latency, per-job cost estimate |
| Spend | Provider dashboards (Google AI Studio usage page, fal.ai billing) — bookmark both |

**Cannot be verified from a phone (accept and defer):** multi-device browser matrix (you have
one physical device), Lighthouse/perf profiling at depth, screen-reader a11y passes, load
testing, and any native-app work. None of these block a demo.

## 3. Ordered task list

Tasks are sized for one agent session each, with a machine-checkable DoD. `[P]` = can run in
parallel with other `[P]` tasks in its phase. Unmarked = serial.

### Phase A — Foundation (days 1–4, serial spine)

- **A1. Scaffold**: Next.js 15 + TS + Tailwind; deploy to Vercel; PWA manifest + icons +
  installability. DoD: preview URL renders "hello" and passes CI; installable from Safari
  share sheet.
- **A2. Supabase project + schema v1**: `profiles`, `garments`, `avatars`, `jobs`, `outfits`,
  `products`, `consents` + RLS on everything + Storage buckets (`uploads`, `garments`,
  `renders`) with owner-only policies. DoD: migration applies in CI; RLS verified by a test
  that a second user cannot read user 1's rows.
- **A3. Auth**: Supabase email OTP (no passwords — phone-friendly). DoD: sign in from preview
  URL on phone.
- **A4. CI + screenshot loop**: Playwright against preview URL, screenshots committed to PR.
  DoD: PR shows fresh screenshots of `/`, `/closet`.
- **A5. Health/selftest route + keepalive workflow.** DoD: `/api/health` green on phone.

### Phase B — Two pipelines (days 5–14, two parallel tracks after A)

Track 1 — Avatar & try-on:
- **B1. Upload flow**: camera/photo-library capture → client-side downscale (~1536px) → signed
  upload. DoD: file lands in `uploads` bucket from phone camera.
- **B2. Provider layer**: `lib/providers` with Gemini image (default) + fal adapter; `jobs`
  create/poll routes; cost-per-job estimate written to the row. DoD: `/api/health` runs a
  test generation end-to-end.
- **B3. Canonical avatar builder**: 1–3 user photos → standardized full-body studio portrait
  (house style prompt v1, fixed background/lighting/crop) → user picks 1 of N → stored as
  canonical. DoD: avatar visible in app for seeded user.
- **B4. Dressing room v1**: pick 1 garment from closet → single-edit try-on from canonical
  avatar → result viewer with save-to-outfits. DoD: end-to-end try-on on preview URL in
  <30s with progress UI.
- **B5. [P] Multi-garment outfits** (one multi-image call) + regenerate button + generation
  history on the job row.

Track 2 — Closet:
- **B6. [P] Garment capture**: photo → client-side background removal → cutout preview →
  accept/retake. DoD: cutout PNG in `garments` bucket, works on phone camera photos of
  clothes on a bed/hanger.
- **B7. [P] Attribute extraction**: Gemini Flash with a strict JSON schema (category,
  subcategory, colors, pattern, fabric guess, formality, season, brand-if-visible) + editable
  form (human-in-the-loop correction is a feature, not a fallback). DoD: 10 test garments
  auto-tagged; edits persist.
- **B8. [P] Closet grid UI**: filterable by category/color; empty states; delete (cascades to
  storage). DoD: screenshots in PR.

### Phase C — The demo layer (days 15–22, mostly parallel)

- **C1. [P] Art direction system**: 2–3 style presets as versioned prompt templates
  (`lib/style/presets.ts`), fixed 4:5 crop, film-grain + grade overlay in CSS, consistent
  typography. Regenerate = same seed policy documented. DoD: same garment rendered in all
  presets, side-by-side page for comparison.
- **C2. [P] Lookbook pre-generation**: `workflow_dispatch` batch job that renders the seeded
  demo closet into ~30 curated outfit images (generate N, keep best via founder review page).
  DoD: `/lookbook` gallery loads instantly (pre-generated, no live inference).
- **C3. [P] Recommendations v1**: rules + Gemini text — "complete this outfit" from closet
  attributes; occasion picker; optional weather (Open-Meteo). No ML training. DoD: 3 sensible
  suggestions for seeded closet, each renderable via try-on.
- **C4. [P] Shop rail**: `products` table from a curated JSON (later: affiliate feed import
  when a network approves); "similar to your closet" via attribute match; outbound affiliate
  links with subid tracking + click log table. DoD: clicking through from phone lands on
  retailer with tag present.
- **C5. [P] Consent & privacy plumbing**: consent screen at upload (explicit biometric-style
  language), `consents` row with timestamp + text version, self-serve "delete everything"
  that cascades DB + Storage, `/privacy` page. DoD: delete leaves zero rows/objects (verified
  by selftest).
- **C6. Latency masking & polish** (after B4): skeleton shimmer, progress states, graceful
  provider-failure retry, PWA offline shell for closet browsing.

### Phase D — Demo hardening (days 23–30, serial)

- **D1. Seeded demo account** with polished closet (25–40 garments), canonical avatar, filled
  lookbook — one-tap "view demo" login.
- **D2. Landing page** with the 5 best renders + install instructions + waitlist email capture.
- **D3. Failure-mode sweep**: airplane-mode test, mid-generation refresh, oversized upload,
  duplicate taps.
- **D4. Demo script**: 90-second walkthrough you can screen-record on the phone; a
  `/press` page with the best stills.
- **D5. Buffer.** Something in B/C will slip. This is that time.

### Parallelization rules (what breaks if run concurrently)

- **Safe in parallel:** Track 1 vs Track 2 in Phase B; all `[P]` tasks in C (they touch
  disjoint routes/dirs).
- **Serial, will break if parallelized:**
  - **Migrations** — one agent owns `supabase/migrations/` at a time; parallel schema edits =
    merge conflicts in ordered SQL files and a broken shadow DB. Queue schema changes.
  - **Provider layer (B2)** must merge before B3/B4/C2 start (they all consume its interface).
  - **Style presets (C1)** before lookbook batch (C2) uses them.
  - Anything editing `app/layout.tsx`, global CSS tokens, or auth middleware — single-writer.
- Practical rule from a phone: run at most 2–3 agent tracks at once, each on its own branch
  and directory; merge trains through PR previews, not simultaneous pushes to one branch.

### Week-by-week milestones

| Week | Milestone — each ends with something you can tap on a phone |
|---|---|
| 1 | Deployed installable PWA shell, auth, schema, CI screenshot loop, health check green. First test generation via API works. |
| 2 | You (founder) have a canonical avatar; 10+ real garments cataloged with cutouts + attributes; first end-to-end try-on on the preview URL. |
| 3 | Art-directed results in 2–3 presets; pre-generated lookbook; recommendations + shop rail with live outbound links; consent/delete flows. |
| 4 | **Demo-ready state:** stranger can open URL → view seeded demo (instant lookbook) → or sign up, build avatar, add 3 garments, run 2 live try-ons — reliably, in under 5 minutes, on mobile Safari. Landing page + 90s recorded demo. |

**Definition of "demo-ready":** the live-generation path works ≥9/10 attempts at <30s; the
wow-path (lookbook) is pre-generated and instant; nothing in the demo depends on a laptop.
