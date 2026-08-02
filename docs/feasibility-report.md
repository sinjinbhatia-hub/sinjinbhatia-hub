# Wardrobe Try-On App — Technical Due Diligence & 30-Day Build Plan

**Date:** 2026-08-02 · **Constraints:** $500 once, 30 days, phone-only via Claude Code on a remote repo, mobile-web PWA, solo founder (Python/FastAPI, React, Supabase, Vercel, Railway known).
**Companion file:** [`BUILD_PLAN.md`](./BUILD_PLAN.md) — the agent-executable task list (Deliverable 4 in full).

**How to read confidence tags:** claims were verified 2026-08-02 by web research. ✅ = corroborated by official docs/pages or multiple independent sources. ⚠️ = secondary sources, conflicting reports, or working-from-memory — re-verify before depending on it. Free-tier quotas are the least stable numbers in this document; re-check them the day you create accounts.

---

## Executive summary — the 5 things that determine whether this ships

1. **The avatar is a 2D editorial illusion, and the identity problem is an architecture decision, not a model purchase.** The only viable pipeline at your constraints: generate ONE canonical avatar image per user (curated, in a fixed house style), then apply garments to it with a dedicated try-on model (FASHN/Kling-class) that inpaints only the garment region and leaves the face pixels untouched. Never chain edits. Get this wrong and every generation is a slightly different person; get it right and identity is solved *by construction* for the core flow.

2. **Money is not the constraint — production discipline is.** Try-on images cost $0.04–$0.075 each; $500 buys ~6,600–12,800 of them, and free tiers (Gemini free tier + Modal's $30/mo credits + a likely $300 GCP trial) can cover an entire demo. The difference between "AI slop" and "cinematic" is curation ratio, one consistent visual treatment, and lighting/crop discipline — all free, all time.

3. **Legal is the real gate to a public link, not tech.** Face photos put you in BIPA/CUBI/My-Health-My-Data territory (private rights of action, $1k–$5k statutory damages per person in IL; Texas AG actively hunting; WA covers face imagery even for non-health apps). A private demo with yourself + consenting adult friends is fine *now* with an e-signed release and a retention policy. A stranger-facing link requires: LLC, ToS, privacy + biometric consent + retention/deletion policy, 18+ gate, and (cheapest risk control) geo-gating IL/TX/WA. Also: **never run try-on against retailer product images** — Amazon's agreement explicitly bans modifying content *and* using it for ML, and no network's feed terms grant derivative-image rights.

4. **The window is real but closing in quarters, not years.** Google killed its standalone try-on app (Doppl, dead April 2026) — try-on survives as a free feature of Google Search. Doji ($14M, Thrive) is still invite-only with accuracy complaints. Whering (10M users) raised $7M in **July 2026** — from eBay and *Google's AI Futures Fund* — explicitly to add try-on to a closet app. The full loop (closet + own-likeness try-on + recs + affiliate checkout) is unshipped by anyone, but the wedge is being funded at someone else right now. Your only edges are speed, taste, and a niche.

5. **Phone-only agentic build is genuinely feasible** — Vercel preview URLs, Supabase dashboard, GitHub Actions `workflow_dispatch` buttons, and CI screenshot artifacts give you a complete verify-from-Safari loop. Risk concentrates in the week-2 generation pipeline and in free-tier instability, not in infra. The stack that survives the constraints is Next.js on Vercel + Supabase — with two license landmines avoided (see Kill List #6 and the Vercel non-commercial clause).

**The honest conclusion you asked for, stated plainly:** the photorealistic 3D avatar is out of reach — not just for you, for everyone at consumer cost (evidence below, incl. Ready Player Me's shutdown and the fact that Doji itself is 2D diffusion wearing 3D branding). Build the 2D diffusion version with strong art direction. That is not the consolation prize; it is what the funded competitors ship.

---

## Deliverable 1 — The $0 demo

### Approach (exact)

**Canonical-avatar + garment-inpainting, two-model pipeline:**

1. **Avatar creation (one-time per user):** 1–3 phone photos → Gemini image model generates a standardized full-body portrait in your house style (fixed pose vocabulary, background, lighting, 4:5 crop). Generate 4–8 candidates; the user picks one. That image is the **canonical avatar** — the single source of truth for their likeness.
2. **Try-on (every generation):** canonical avatar + garment cutout → try-on model → result. Dedicated try-on models (FASHN v1.6, Kling Kolors) work by inpainting the garment region onto your input photo — the face is *copied, not regenerated*, so likeness holds by construction ✅ ([FASHN docs](https://docs.fashn.ai/api-reference/tryon-v1-6), [fal Kling page](https://fal.ai/models/fal-ai/kling/v1-5/kolors-virtual-try-on)). General editors (Gemini) regenerate the whole image — flexible, but identity drifts.
3. **Rule that must never be broken:** every try-on is a single hop from the canonical image. Never edit an edited image — drift compounds ✅ (practitioner consensus; see D3).
4. **Art direction layer:** versioned prompt templates + fixed crop + CSS grain/grade overlay. Free. This is where "cinematic" comes from (D3).

At $0 specifically, you can't use FASHN (no free tier), so the $0 demo runs on Gemini free tier for both steps and accepts curation overhead — with Modal free credits as a batch fallback. This is the first honest tension: **the $0 engine is not the engine you'd choose with $40.** See D2.

### Free-tier inference options, verified 2026-08-02

| Provider / model | Per-image (paid) | Free allowance | Notes |
|---|---|---|---|
| **Gemini 2.5 Flash Image** ("Nano Banana") | $0.039 (≈1,290 output tokens @ $30/M); ~$0.0195 batch ⚠️ | Free tier still exists; reports conflict: up to ~500 images/day claimed vs. post-Dec-2025 quota cuts → **plan on ~100/day** ⚠️ | Default engine at $0. ⚠️ One source reports retirement Oct 2026 — check migration path. Free-tier data may be used by Google to improve models ⚠️ (long-standing terms; verify at signup — matters for friends' photos). |
| **Gemini 3 Pro Image** ("Nano Banana Pro") | $0.134 (1–2K), $0.24 (4K); batch ~50% off ⚠️ | **No free API tier** ⚠️ | Hero-shot tier. On fal: $0.15/img ⚠️. |
| **fal.ai — FASHN try-on v1.6** | $0.075 ✅ | No meaningful free credits (signup credits exist, amount unverified ⚠️) | Best garment-fidelity API; 864×1296 native; 5–17s ✅. |
| **fal.ai — Kling Kolors try-on v1.5** | $0.07 ✅ | — | Cheaper alternative; holds pose/body/skin ✅. |
| **fal.ai — BiRefNet v2 (bg removal)** | $0.0008/compute-sec (≈ sub-cent/image) ✅ | — | For cataloging (D5). |
| **Replicate — idm-vton** | ~$0.025/run ⚠️ | **No free signup credits** (prepaid model since Jul 2025) ✅ | Model is CC BY-NC — non-commercial (D3). |
| **Modal** | GPU: T4 ~$0.59/hr, A10G ~$1.10/hr, A100 ~$2.10–3.40/hr ⚠️ | **$30/month free compute credits** on the $0 Starter plan ✅ | Run **Leffa** (MIT license, ~6s/img on A100) as a batch try-on engine: ≈$0.01–0.02/img amortized → **~1,500–3,000 free images/mo**. Batch only (cold starts); Python — your language. |
| **Hugging Face ZeroGPU** | — | ~3.5 GPU-min/day free; PRO $9/mo → ~25 min/day ✅ | Dev toy / experimentation, not demo infra. |
| **Google AI Studio (web UI)** | — | Interactive use free | Manual pre-generation of lookbook assets costs $0. |
| **Vertex AI `virtual-try-on-001`** | GA, price not publicly verifiable (~$0.04 est. ⚠️) | **New GCP accounts historically get $300/90-day trial credits** ⚠️ verify | If the $300 trial still exists, it alone funds the whole demo on paid-tier privacy terms. Check first. |

Sources: [fal FASHN](https://fal.ai/models/fal-ai/fashn/tryon/v1.5) · [fal Kling](https://fal.ai/models/fal-ai/kling/v1-5/kolors-virtual-try-on) · [fal BiRefNet](https://fal.ai/models/fal-ai/birefnet/v2) · [Gemini pricing tracker](https://pricepertoken.com/pricing-page/model/google-gemini-2.5-flash-image) · [free-tier reports](https://www.aifreeapi.com/en/posts/gemini-api-free-tier-rate-limits) · [NB Pro pricing](https://www.aifreeapi.com/en/posts/gemini-3-pro-image-preview-pricing-calculator) · [Replicate free-plan status](https://costbench.com/software/ai-productivity/replicate/free-plan/) · [Modal pricing](https://www.beam.cloud/blog/modal-pricing-explained) · [HF ZeroGPU](https://huggingface.co/docs/hub/en/spaces-zerogpu) · [Leffa (MIT)](https://github.com/franciszzj/Leffa) · [Vertex try-on GA](https://docs.cloud.google.com/vertex-ai/generative-ai/docs/model-reference/virtual-try-on-api)

**Notable churn found during research:** fal appears to have removed its idm-vton endpoint; Imagen 4 is deprecated (shutdown Aug 17, 2026); Replicate ended free credits; a "Gemini 3.1 Flash Image / Nano Banana 2" preview (~$0.05/img) is reported but weakly sourced ⚠️.

### How many generations does $0 actually buy?

- **Live:** ~100/day on Gemini free tier (conservative planning number; could be up to 500/day ⚠️) → ~3,000/month.
- **Batch:** ~1,500–3,000/month via Leffa on Modal's $30/mo free credits (setup cost: ~1–2 agent-days; it's Python).
- **Manual:** unlimited-ish interactive generations in AI Studio's web UI for pre-built lookbook assets.
- **Practical answer: free capacity is not the binding constraint for a 1–10 user demo.** The binding constraints are identity consistency and curation time.

### What the $0 demo can and cannot show a viewer

**Can show:** an installable PWA opened from a link; a real cataloged closet (cutouts + attributes); a canonical avatar that looks like the user in a consistent editorial style; async try-on of owned garments in ~10–25s with progress UI; an instant, pre-generated 20–30-image lookbook (the wow moment — because it's curated); outfit recommendations; a shop rail with real product links (images unmodified).

**Cannot show:** guaranteed likeness on every live generation (Gemini free tier regenerates the whole image — expect visible identity variance strangers will notice); garment-critical fidelity (logos, prints) at FASHN level; video; real-time anything; try-on of arbitrary shoppable items (legal, D6/D7); revenue.

---

## Deliverable 2 — What $500 changes

### Recommended allocation

| Line item | $ | What it buys | Marginal return |
|---|---|---|---|
| FASHN v1.6 / Kling budget (via fal) | $50 | ~600–700 garment-faithful try-ons | **Highest.** Swaps the demo's core engine from "whole-image regen with drift" to "face-preserving inpainting." First $50 changes the demo more than the last $450. |
| Gemini paid tier (avatar gen + cataloging + recs) | $100 | ~2,500 images at $0.039 + effectively-free attribute extraction | **High.** Volume for curation (8:1 generate-to-keep) + removes free-tier rate anxiety + paid-tier data terms (photos not used for training ⚠️ verify) — matters the moment you process friends' faces. |
| Nano Banana Pro hero shots | $30 | ~220 images at $0.134 | **Medium.** Landing page / lookbook heroes only. Identity degrades on full-body garment swaps (research below) — use for styled portraits, not try-on. |
| Domain (.com, Porkbun/Cloudflare) | $12 | Credibility | High per dollar; trivial. |
| Awin deposit | $5 | Network access + product feeds (refundable) | High per dollar (D6). |
| Founder-only identity LoRA experiment | $10 | 3–5 FLUX LoRA training runs ($2/run fal, ~$1.50 Replicate ✅) | **Low-medium.** Worth one experiment for pitch assets; not the product architecture (see below). |
| Vercel Pro, 1 month — when affiliate links go live | $20 | Hobby tier prohibits commercial use ✅ | Required *at monetization*, not before. (Alternative: move to Cloudflare Pages, $0.) |
| LLC reserve (home state) | $150 | The public-launch legal gate (D7) | Zero demo value; unlocks a public link. Spend only when going public. |
| Buffer / overruns | $123 | Reality | — |
| **Total** | **$500** | | |

### The comparisons you asked for

- **More inference volume:** $500 all-in buys ~12,800 Gemini images ($0.039), ~7,100 Kling ($0.07), ~6,600 FASHN ($0.075), ~3,700 NB Pro ($0.134). You will not exhaust this in 30 days of solo building — volume is the *wrong* thing to max.
- **Better model tier:** NB Pro is ~3.4× Gemini Flash cost and *worse* at full-body garment replacement (identity degradation, per an Apr 2026 arXiv comparison ⚠️ competing-vendor paper). Model tier upgrades buy polish on hero shots, not a better product. **Cap at ~$30.**
- **Fine-tune / LoRA of one person:** $2–4 and ~20 min per person, real identity gains for *generative* shots — but it doesn't solve garment fidelity, doesn't scale operationally (a training run per user), and the canonical-avatar + inpainting pipeline gets you most of the benefit for free. **Do one for yourself for the pitch; don't build the product on it.**
- **Freelancer for 3D garment assets or art direction:** $500 buys roughly 3–10 game-grade 3D garments (useless without a 3D pipeline, which is dead — D3) or ~5–15 hrs of junior art direction. Art direction is the one skill worth *learning* instead — it's prompt templates + reference boards, and agents can execute a written style guide. **$0.**
- **Domain / Apple dev / entity:** Domain yes ($12). Apple Developer ($99) — no: you have no Mac, no Xcode, and a PWA is the deliberate architecture; buying it is pure waste this month. Entity — yes as a *reserve*, triggered by going public, not by the demo.

**Verdict:** $500 makes the demo meaningfully better in exactly two ways — the FASHN/Kling engine swap (identity + garment fidelity) and paid-tier privacy for other people's photos. Roughly $200 of it is genuinely productive this month; the rest is a launch reserve. If someone forced you to build at $0, you'd still have a credible demo; if they gave you $5,000, the avatar would look barely different — **spend threshold ~$200, save the rest.**

---

## Deliverable 3 — The quality ceiling (unsoftened)

### What the avatar actually looks like at $500 / 30 days

Static stills, 864×1296 to 2K, 4:5 crop. On the FASHN/Kling pipeline: **the face is literally your photo** (inpainting preserves it), the garment is ~85–95% faithful (logos and fine print mostly right at v1.6's resolution, occasionally mangled), and the tells are fabric drape physics, hands, hair-garment boundaries, and lighting coherence between body and background. On the Gemini pipeline (avatar creation, stylized shots): "a person who strongly resembles you" — the sibling problem — right ~80–90% of single-hop generations, with occasional off-model faces you must curate out. Latency 5–25s. No video (except 1–3 pre-rendered garnish clips). Body shape is approximate and tends to idealize — a product-trust and ethics issue you must decide on deliberately (Doji gets rated 3.4★ partly for making people thinner/taller ⚠️ [App Store reviews via comparison sites](https://www.drapeai.app/drapeai-vs-doji)).

**The research reality check:** an April 2026 arXiv benchmark reports Nano Banana Pro, GPT-Image-1.5, and Seedream-5-lite all suffer "severe identity degradation" on complex full-body garment replacement vs. dedicated try-on models ⚠️ ([arXiv 2604.19748](https://arxiv.org/pdf/2604.19748) — a competing model's paper, but consistent with practitioner reports). This is why the two-model architecture in D1 is the load-bearing decision of the entire project.

### The gap between that and "ultra-realistic, cinematic"

What that phrase is, in every shipped example you've seen (Doji's marketing, fashion-AI campaign work):

1. **Curation ratio** — generate 8–20, publish 1. The reel is the top 5%. Cost: time, ~$0.50/kept image.
2. **One consistent treatment** — a single lighting scheme, palette, crop, film-grain pass, and typeface across every image. Cost: $0, one day of writing a style guide. This is the single biggest perceived-quality lever you control.
3. **Model tier on heroes only** — NB Pro/4K for the landing page. Cost: ~$30.
4. **Motion garnish** — image-to-video on 1–3 best stills. Cost: a few dollars.
5. **Human retouch** — what agencies add. Skip.

**Answer to "is art direction cheaper than realism?" — Yes, by orders of magnitude, and it's also the only one of the two you can actually buy.** Realism's asymptote is set by the frontier models you rent; $500 and $50,000 rent the same models. Art direction is taste + consistency + curation, all available at $0. A stylized, confident "house look" also *hides* AI artifacts that a neutral photoreal look exposes — leaning artsy is the strictly dominant strategy at your budget.

### The 3D avatar path — proven dead (as requested)

- **Ready Player Me: gone.** Acquired by Netflix Dec 2025; public platform and APIs shut down Jan 31, 2026 ✅ ([Variety](https://variety.com/2025/digital/news/netflix-acquires-ready-player-me-games-avatar-creation-1236612915/)). The "free interoperable avatar" path no longer exists.
- **Avaturn:** realistic *head* on parametric game-grade bodies; pricing floor ~$800/mo ⚠️ ([pricing page via search](https://avaturn.me/pricing/)). Not photoreal, not in budget.
- **Meshcapade / SMPL-X:** SMPL body models are research-only; commercial licensing is a negotiated enterprise contract with Meshcapade — no public price ✅ ([license](https://smpl-x.is.tue.mpg.de/modellicense.html), [Meshcapade](https://meshcapade.com/smpl/)). Produces accurate *bodies*, not photoreal appearance.
- **Gaussian splatting:** no hosted service anywhere turns a phone capture into a photoreal *riggable* avatar in 2026 ✅ (research pass found none); pipelines are papers requiring per-subject GPU optimization and careful capture ([SplattingAvatar](https://ieeexplore.ieee.org/iel8/10654794/10654797/10656303.pdf), [arXiv 2510.13978](https://arxiv.org/pdf/2510.13978)). And even a perfect splat of you doesn't solve try-on — *re-dressing* a splat is an open research problem, and you'd still need rigged 3D garment assets that don't exist for your closet.
- **The clincher:** Doji — $14M from Thrive, founders ex-Apple VisionOS and ex-DeepMind — ships **2D diffusion under 3D-avatar branding** ✅ ([TechCrunch](https://techcrunch.com/2025/05/15/doji-raises-14m-to-make-virtual-try-ons-fun-through-ai-avatars)). If photoreal 3D try-on were shippable, the team best equipped on Earth to ship it would not be doing image generation.

**Conclusion: 2D is not your budget compromise; it is the industry's current ceiling for this product.** Your disadvantage vs. Doji is their in-house model tuning and $14M of curation/eval labor — which you partially offset by renting FASHN (garment fidelity is *their* whole business) and out-taste-ing them in one niche.

---

## Deliverable 4 — 30-day build plan for agentic execution

**The full ordered task list, DoD per task, and parallelization rules live in [`BUILD_PLAN.md`](./BUILD_PLAN.md).** Summary of the decisions and why:

- **Stack:** Next.js 15 on Vercel Hobby + Supabase Free + GitHub Actions. One deployable, preview URL per PR (your phone verification loop), API routes hold secrets, and it's the stack coding agents execute most reliably — which matters more than your FastAPI familiarity when you can't debug locally. Python survives in exactly one place: an optional Modal batch worker (Leffa) deployed from CI.
- **Free-tier landmines designed around (verified ✅):** Vercel Hobby functions cap at 300s (Fluid) — fine, because generation is client-polled jobs, never long-held requests; Vercel Hobby is **non-commercial only** — fine for the demo, swap to Pro ($20/mo) or Cloudflare Pages when affiliate links go live; Supabase Free pauses after ~7 idle days — a scheduled GitHub Actions ping prevents it; Supabase Edge Functions 150s / uploads 50MB / storage 1GB — all sufficient at demo scale; GitHub Actions 2,000 min/mo private (make the repo public → unlimited).
- **Verification from a phone:** Vercel preview URLs; Playwright screenshots committed to the PR (GitHub mobile renders them); Supabase dashboard; a `/api/health` selftest route; `workflow_dispatch` buttons in the GitHub mobile app for seeding, lookbook batch generation, and smoke tests; provider spend dashboards bookmarked.
- **Cannot verify from a phone (defer, doesn't block demo):** multi-device browser matrix, deep performance profiling, screen-reader passes, load testing.
- **Parallelization:** after the week-1 spine (scaffold → schema → auth → provider layer), Track 1 (avatar/try-on) and Track 2 (closet cataloging) run in parallel, then four parallel week-3 workstreams (art direction, lookbook batch, recommendations+shop, consent/privacy). **Serial and will break if parallelized:** Supabase migrations (single writer), the provider abstraction (everything consumes it), global layout/auth middleware. Practical cap from a phone: 2–3 agent tracks at once, merged through PR previews.
- **Milestones:** W1 installable PWA + auth + CI loop + first API generation; W2 founder's canonical avatar + 10 cataloged garments + first end-to-end try-on; W3 art-directed presets + pre-generated lookbook + recs + shop rail + consent/delete flows; W4 hardening + seeded demo account + landing page + 90s screen-recorded demo. **Demo-ready definition:** live path works ≥9/10 at <30s; wow path is pre-generated and instant; nothing requires a laptop.

---

## Deliverable 5 — Closet cataloging on a budget

Pipeline per garment: phone photo → background removal → cutout → attribute extraction → editable form → catalog row.

| Step | Recommended | Cost | Alternatives & license traps |
|---|---|---|---|
| Background removal | **fal BiRefNet v2** server-side | $0.0008/compute-sec → **~$0.001–0.003/image** ✅ | **@imgly/background-removal is AGPL-3.0** — poison for a closed commercial app ✅ ([LICENSE](https://github.com/imgly/background-removal-js/blob/main/LICENSE.md)). **BriaAI RMBG-1.4/2.0: non-commercial** ✅. License-clean $0 option: **BiRefNet weights are MIT** ✅ ([repo](https://github.com/ZhengPeng7/BiRefNet)) — run in-browser via transformers.js/ONNX (WebGPU shipped in iOS 26 ⚠️), at the cost of ~1–2 agent-days and a 40–90MB model download. Start with fal (pennies), add in-browser later. |
| Attribute extraction | **Gemini 2.5 Flash** (vision, strict JSON schema): category, subcategory, colors, pattern, fabric guess, formality, season, brand-if-visible | ~$0.0001–0.0005/image paid; free tier covers demo ✅ | Fashion-CLIP embeddings later for similarity; not needed for v1. |
| Correction UI | Editable form, always | $0 | Non-negotiable: treat auto-tags as drafts. |

**Total: well under half a cent per garment.** A 40-item closet costs ~$0.10 to catalog. Accuracy expectations (be honest in the UI): top-level category ~90–95%, dominant color high-90s, pattern good, **fabric/material ~70%** (visual guessing), brand only when a logo/label is visible. The human-in-the-loop edit step is what makes the catalog trustworthy — and it doubles as engagement.

Photography guidance in-app (flat on a bed, hanger on a door, decent light) moves accuracy more than any model choice.

---

## Deliverable 6 — Commerce, realistically

### Who approves a zero-traffic pre-launch app (verified 2026-08-02)

| Network | Zero-traffic approval? | Cost | Product feed w/ images? | Notes |
|---|---|---|---|---|
| **Sovrn Commerce** | Yes — lowest bar, days ✅ | $0, 75/25 split | Product/price APIs ✅ (tiny-account tier ⚠️) | **Start here.** H&M US reachable ⚠️. |
| **Awin** | Yes — accessible to new publishers ✅ | $5 refundable deposit ✅ | **Best**: Create-a-Feed + "soft membership" feeds before per-advertiser approval ✅ | ASOS (EU/UK) here ⚠️. **ShareASale is gone** — sunset Oct 6, 2025, merged into Awin ✅. |
| **Impact.com** | Marketplace review ~2 days; some brands auto-approve ✅ | $0 | Catalogs via UI/FTP/API ✅ | Third to join. |
| **Rakuten** | Wants a live site with content + traffic numbers ✅ | $0 | Per-advertiser | After launch. Nordstrom, SSENSE, ASOS-US ⚠️. |
| **CJ** | Hardest — per-advertiser gates (site age, content volume) ✅ | $0 | Per-advertiser | Later. Revolve ✅, Uniqlo ⚠️. |
| **Skimlinks** | Operating (Taboola) but **2025–26 reports of withheld payments** ⚠️ | 75/25, $65 min payout | API historically | Avoid for now. |
| **Amazon Associates** | Provisional instant, **but**: 3 sales in 180 days or auto-withdrawn ✅; **PA-API retired May 2026**, successor Creators API needs **10 qualified sales / trailing 30 days** ✅; **Operating Agreement bans altering images and using content for ML** ✅ | $0, apparel 4%, 24h cookie | Effectively no (pre-traction) | **Wrong first network in 2026.** ([Creators API gate](https://blog.freshstore.com/amazon-creators-api-pa-api-retirement/), [OA](https://affiliate-program.amazon.com/help/operating/agreement)) |

### Feeds, free catalogs, and the try-on trap

- **Free/open fashion catalogs are all non-commercial:** H&M Kaggle dataset, DeepFashion/2 — research-only ✅; Fashionpedia — unclear/Flickr-sourced ⚠️. MIT-licensed **DummyJSON/FakeStoreAPI** work for UI demos but look generic ✅. There is **no commercial-safe open fashion catalog.** Cleanest real-catalog path: Awin $5 → soft-membership feeds.
- **The trap:** no network's standard feed terms grant the right to *modify* product images, and Amazon's explicitly prohibits it (plus an ML-use ban) ✅. **Do not run try-on against feed/retailer images.** The compliant pattern: try-on runs on the user's own garments; the shop rail shows retailer images **unmodified** with "similar to your closet" matching. Try-on of *purchasable* items needs a written license or a direct brand partnership — that's a post-demo BD task, and honestly a wedge: brands might pay for it.

### Commission math (why revenue is a non-goal this month)

Fashion pays 4–10% (Amazon 4%, Revolve 5%, ASOS ≤7%) with 1–3% conversion on referral clicks, $80–120 AOV, and 10–25% apparel return/reversal rates ✅⚠️. **Per 1,000 product clicks: ~$35–$230 net (EPC ≈ $0.04–0.23).** A demo with 50 users generating a few hundred clicks earns approximately nothing. Build the click-tracking plumbing (subid per user/outfit) because it's cheap and proves the loop; treat actual revenue as a 2027 problem.

---

## Deliverable 7 — Legal gates before any public launch

*(Fact-collection with citations, not legal advice. The one paid professional worth consulting before a public link is a privacy lawyer for one hour.)*

### What you can do NOW (private demo, yourself + consenting friends)

Legally lowest-risk configuration, achievable this week: adults only (18+), each tester signs an **electronic consent/release** before upload (BIPA's 2024 amendment SB 2979 expressly validates e-signatures ✅ [Greenberg Traurig](https://www.gtlaw.com/en/insights/2024/8/bipa-update-illinois-limits-liability-and-clarifies-electronic-consent-for-biometric-data-collection)), a **published retention/destruction policy** (even a markdown page), a working delete-everything button, invite-gated access, and — if any tester is in IL/TX/WA — treat their data with the full consent flow or skip those testers. Send friends' photos through **paid-tier** APIs, not free tier (free-tier data may be used for model improvement ⚠️ verify terms) — this costs dollars.

### What changes the moment a stranger uploads a photo

You become a biometric data controller at scale-zero with class-action surface:

- **Illinois BIPA:** private right of action, $1,000/$5,000 statutory damages per person (SB 2979 killed per-scan stacking, not per-person exposure) ✅. Courts have let "face geometry from photos" claims proceed (*Rivera v. Google*, *Patel v. Facebook* — $650M settlement) ✅; **directly analogous apps have been sued** (Lensa avatar app, Feb 2023; a wave of try-on retailer suits) ✅ ([TFL](https://www.thefashionlaw.com/ai-avatar-generator-lensa-lands-on-receiving-end-of-biometric-data-lawsuit/), [ArentFox Schiff](https://www.afslaw.com/perspectives/alerts/biometric-privacy-class-actions-take-aim-virtual-try-retailers)). Your best factual defense — the pipeline copies face pixels rather than extracting an identifying template (*Martell v. X Corp.* line) — is untested for diffusion try-on ⚠️. Also §15(c): a flat no-profiting-from-biometrics rule with no consent cure — keep monetization attached to recommendations/links, never to the face data itself.
- **Washington My Health My Data:** sleeper risk — "consumer health data" includes biometric data, and commentators read face **imagery** as qualifying even in non-health apps; private right of action via the CPA ✅⚠️ ([Hintze Law](https://hintzelaw.com/blog/wa-my-health-my-data-act-pt7-biometrics)).
- **Texas CUBI:** AG-only but $25k/violation, and the AG's office is the most aggressive biometric enforcer in the country (Meta $1.4B, Google $1.375B settlements) ✅.
- **Colorado (HB 24-1130, eff. July 2025):** biometric consent + written retention/deletion/breach policy **with no user-count threshold** — it applies to you at user #1 ✅ ([Covington](https://www.insideprivacy.com/data-privacy/colorado-privacy-act-amended-to-include-biometric-data-provisions/)).
- **Tennessee ELVIS Act:** liability if users upload *someone else's* photo — you need an "only yourself" attestation + report/remove flow ✅.
- **GDPR (if you serve the EU):** likely *not* Art. 9 special-category as long as you never compute identification/matching embeddings (Recital 51: photos are biometric only when processed for unique identification ✅) — but still personal data (lawful basis, DPIA likely). **EU AI Act Art. 50 applies from Aug 2, 2026 — today:** AI-generated images must carry machine-readable marking/disclosure ⚠️. Add "AI-generated" labeling + metadata regardless; it's cheap and inevitable.
- **Minors:** COPPA's 2025 amendments explicitly add biometric identifiers ✅. 18+ gate, delete-on-sight for suspected minor photos.

### Minimum consent + retention + deletion (the checklist)

1. Separate biometric-consent screen (not buried in ToS): what's collected, why, that AI processes it, which vendors see it (Google/fal), how long it's kept; affirmative e-sign action; stored consent record with policy version.
2. Public retention/destruction policy: originals deleted within N days of avatar creation (pick 7–30); renders user-deletable; hard outer bounds — CUBI 1yr post-purpose, BIPA 3yr from last interaction ✅.
3. Self-serve delete-all that cascades DB + storage + a documented vendor-deletion step; deletion verified by an automated test.
4. No sale/share/training-use of face data, stated and honored; paid API tiers only for user photos.
5. 18+ gate; "photos of yourself only" attestation.
6. **Geo-gate IL, TX, WA at public launch** (the standard cheap startup mitigation ⚠️ common practice) until a lawyer reviews your consent flow — then open them.

### Brand/IP exposure (try-on rendering)

Rendering a real brand's garment **from the user's own closet photo** is low-risk: first-sale/nominative principles cover identifying genuine goods; use brand names in text, no logos in your UI chrome, and a visible "AI simulation — not brand-endorsed" disclaimer ✅ ([INTA](https://www.inta.org/fact-sheets/fair-use-of-trademarks-intended-for-a-non-legal-audience/)). The real exposure is **copyright in product photography**: generating from retailer images = reproduction/derivative-work claims plus affiliate-ToS breach (Amazon explicitly) ✅ — already killed in D6. Trade dress on rendered garments is a tail risk (a Midjourney trade-dress claim survived dismissal ⚠️) — mitigated by accuracy + disclaimers.

### What must exist before a public link

LLC in your home state (~$50–$150 + $0–125/yr agent; Wyoming only helps if you have no home-state nexus ✅ [fee table](https://www.llcuniversity.com/llc-filing-fees-by-state/)) → ToS + privacy policy + biometric consent + retention policy pages → deletion flow → age gate → AI-content labeling → geo-gates. **Insurance (tech E&O/cyber ~$1.5k–5k/yr, and biometric claims increasingly excluded ⚠️ [Vouch](https://www.vouch.us/blog/startup-insurance-costs)) does not fit in $500** — the honest posture: stay invite-gated and geo-gated until revenue or funding pays for it. That is a conscious, disclosed risk acceptance.

---

## Deliverable 8 — Competitive reality

| Player | Shipped | Raised | The lesson for you |
|---|---|---|---|
| **Doji** | Invite-only iOS avatar try-on (6 selfies + 2 full-body → ~30 min avatar build); social sharing; checkout "planned." Still invite-only 15 months post-seed; ~3.4★ on a tiny ratings base; accuracy complaints ("thinner/taller," "doesn't look like me") ✅⚠️ | $14M seed, May 2025, Thrive + Seven Seven Six ✅ ([TC](https://techcrunch.com/2025/05/15/doji-raises-14m-to-make-virtual-try-ons-fun-through-ai-avatars)) | 2D diffusion behind 3D branding; identity quality is hard even with $14M and elite founders. Your D1 architecture dodges their hardest problem. |
| **Whering** | Closet digitization, outfit builder, styling recs, resale; **no try-on yet — the new raise is earmarked to build it** ✅ | **$7M, July 2026, eBay Ventures + Google AI Futures Fund**; ~$14M total; claims 10M registered users ✅ ([tech.eu](https://tech.eu/2026/07/07/whering-lands-7m-as-digital-wardrobe-platform-reaches-10m-users/)) | Your exact wedge, funded, with distribution — ship before their feature does, and differentiate on look/niche. Also: 5 years to $14M total says closet apps alone don't compound. |
| **Acloset (Looko)** | AI closet + styling chat; 4.5M cumulative / 500K MAU claimed ⚠️ | Early rounds; big numbers unverifiable ⚠️ | Cataloging is commodity; retention is the fight. |
| **Indyx** | Closet + human stylist marketplace ✅ | Pre-seed (undisclosed) ⚠️ | Human-in-the-loop is their answer to AI quality. |
| **Style DNA** | Selfie → color analysis → recs; ~$3M/yr subscription claim ⚠️ | Bootstrapped-ish; MWM-published ⚠️ | Solo-scale fashion apps *can* make real subscription money in a niche. |
| **Google** | **Doppl (standalone try-on app): shut down April 30, 2026** ✅ ([support](https://support.google.com/labs/answer/16537062)). Survivor: own-photo "Try it on" inside Search/Shopping — free, billions of listings, US/AU/CA/JP ✅ | — | The most important strategic datapoint: Google concluded try-on is a **feature of shopping search, not a destination app**. Don't compete on generic try-on; own the closet + identity + taste layer Google has no interest in. |
| **Amazon** | Prime Try Before You Buy killed Jan 2025 ✅; shoes-AR + size recs; no own-photo apparel try-on shipped | — | Even Amazon treats try-on as checkout garnish. |
| **Walmart/Zeekit** | "Be Your Own Model" own-photo try-on since 2022 (~$200M acquisition ⚠️); live but visibly stagnant since 2023 ⚠️ | — | A ~$200M acquisition went quiet: try-on alone doesn't move retail metrics. |
| **New entrants** | Alta ($11M, Menlo — AI stylist w/ closet) ✅; Daydream ($50M — fashion search) ✅; Gensmo (claims $60M "angel" ⚠️); Phia ($43M total — price comparison) ✅ | | Capital is flooding adjacent layers; nobody owns your exact combination. |

**The gap that's actually open:** no funded player ships closet cataloging **+** own-likeness try-on **+** recommendations **+** affiliate checkout as one loop. Whering has three minus try-on (funded to add it); Doji has try-on minus closet; Google has try-on + commerce minus any interest in your existing wardrobe; Alta has closet + recs minus photoreal try-on. **But be clear-eyed about *why* it's open:** the giants concluded try-on is a feature, not a company, and the leader in your wedge just took Google-affiliated money to close it. The defensible version of this for a solo founder is not "the full loop, general market" — it's the full loop for **one niche with one strong aesthetic**, shipped in weeks, as a wedge/portfolio/learning vehicle with a real (small) chance of catching a current.

---

## KILL LIST — not feasible at $500 / 30 days / phone-only

| # | Killed | Why (evidence) | Cheapest credible substitute |
|---|---|---|---|
| 1 | Photoreal 3D avatar (any flavor: RPM, Avaturn, SMPL, splats) | RPM dead ✅; Avaturn ~$800/mo ⚠️; SMPL = enterprise negotiation ✅; splats = papers, and re-dressing them is unsolved ✅; Doji itself is 2D ✅ | Canonical 2D avatar in a fixed house style (D1) |
| 2 | Real-time / interactive dressing | Physics + latency don't exist at any rentable price | Async 5–25s generation + pre-generated instant lookbook |
| 3 | Try-on against retailer/feed/Amazon product images | Copyright/derivative-work exposure + Amazon OA bans modification *and* ML use ✅; no network grants image-modification rights ✅ | Try-on on user's own garments; shop rail shows retailer images unmodified; brand-licensed try-on = future BD |
| 4 | Per-user LoRA identity fine-tuning as product architecture | $2–4 + ~20 min/user, unscalable ops, unpredictable quality ✅ | Canonical-photo single-hop edits + garment inpainting; one founder LoRA ($10) for pitch assets only |
| 5 | Native iOS app | No Mac/Xcode; $99; App Store review time | PWA (installable, push-capable on iOS 16.4+ ✅) |
| 6 | @imgly/background-removal & BriaAI RMBG in the product | AGPL-3.0 and non-commercial licenses respectively ✅ | fal BiRefNet v2 (pennies) or MIT BiRefNet in-browser ✅ |
| 7 | Amazon Associates as the launch network | 3-sales/180d gate + Creators API needs 10 sales/30d ✅ | Sovrn → Awin ($5) → Impact; Amazon post-traction |
| 8 | Public open signup at day 30 | Biometric-law surface with no entity/policies/insurance (D7) | Invite-gated waitlist; geo-gate IL/TX/WA; LLC + policies then public |
| 9 | Video try-on as a core feature | Cost/latency/quality; identity drift worst in video | 1–3 pre-rendered hero clips as marketing garnish |
| 10 | Insurance this month | E&O/cyber ~$1.5–5k/yr, biometric exclusions common ⚠️ | Gate exposure (invites + geo) and accept the residual risk knowingly |
| 11 | Training/fine-tuning your own try-on model | Every open model is non-commercial except Leffa (dataset-provenance risk ⚠️) ✅; no GPU, no ML background, no time | Rent FASHN/Kling; optional Leffa-on-Modal for free batch volume |

---

## Open questions to answer before writing code

1. **Which US state are you (and your first testers) in?** IL/TX/WA/CO residency changes the consent flow and whether geo-gating matters on day 1. (Cheap to answer; changes D7 posture.)
2. **Who is the demo for** — recruiters/investors/friends (stay private, $0 legal) or a public social launch (triggers the full D7 gate + ~$200 spend)? Pick now; it sets the week-4 scope.
3. **One niche + one aesthetic.** "Cinematic" must resolve to a single writable style guide (e.g., quiet-luxury editorial vs. Y2K flash vs. streetwear archive). The art-direction layer needs ONE look — which?
4. **Day-1 free-tier verification:** create the Google AI Studio key and confirm the image model's actual free RPD (conflicting reports: ~100 vs 500/day ⚠️); confirm whether the $300 GCP trial credit still exists ⚠️; confirm fal signup credits ⚠️. These three facts set the real $0 capacity — they're 20 minutes of checking from a phone.
5. **Do you have seed content** — 20–40 of your own garments photographable this week, plus 2–3 consenting adult friends for avatars? The demo's credibility is the seeded closet.
6. **Recommendations v1 scope:** in-closet outfit builder ("wear what you own") or shop-first ("buy this to match")? One is retention, the other is monetization — pick primary.
7. **EU users in scope?** If yes: AI Act Art. 50 labeling (do it anyway) + GDPR posture. If no: block EU at launch and simplify.
8. **Are you comfortable spending ~$10 immediately** (paid Gemini tier) the moment any non-self photo is processed? (Privacy-correct answer: yes.)

---

## Source appendix (primary links)

Pricing/infra: [fal FASHN](https://fal.ai/models/fal-ai/fashn/tryon/v1.5) · [FASHN v1.6](https://fashn.ai/blog/fashn-v1-6-our-best-virtual-try-on-model-yet-now-at-864-x-1296-resolution) · [fal Kling try-on](https://fal.ai/models/fal-ai/kling/v1-5/kolors-virtual-try-on) · [fal BiRefNet v2](https://fal.ai/models/fal-ai/birefnet/v2) · [fal FLUX LoRA training](https://fal.ai/models/fal-ai/flux-lora-fast-training) · [Gemini image pricing](https://pricepertoken.com/pricing-page/model/google-gemini-2.5-flash-image) · [Gemini free-tier changes](https://www.aifreeapi.com/en/posts/gemini-api-free-tier-rate-limits) · [NB Pro pricing](https://www.aifreeapi.com/en/posts/gemini-3-pro-image-preview-pricing-calculator) · [Vertex try-on API](https://docs.cloud.google.com/vertex-ai/generative-ai/docs/model-reference/virtual-try-on-api) · [Replicate prepaid shift](https://costbench.com/software/ai-productivity/replicate/free-plan/) · [Modal pricing](https://www.beam.cloud/blog/modal-pricing-explained) · [HF ZeroGPU](https://huggingface.co/docs/hub/en/spaces-zerogpu) · [Vercel Hobby](https://vercel.com/docs/plans/hobby) · [Vercel image optimization](https://vercel.com/docs/image-optimization/limits-and-pricing) · [Supabase pausing](https://supabase.com/docs/guides/platform/free-project-pausing) · [Supabase function limits](https://supabase.com/docs/guides/functions/limits) · [GitHub Actions billing](https://docs.github.com/en/billing/concepts/product-billing/github-actions) · [Render free](https://render.com/docs/free) · [Railway plans](https://docs.railway.com/pricing/plans) · [Cloudflare R2](https://developers.cloudflare.com/r2/pricing/) · [Open-Meteo terms](https://open-meteo.com/en/terms)

Models/licenses: [IDM-VTON](https://github.com/yisol/IDM-VTON) · [OOTDiffusion license](https://raw.githubusercontent.com/levihsu/OOTDiffusion/main/LICENSE) · [CatVTON](https://github.com/Zheng-Chong/CatVTON) · [FitDiT](https://github.com/BoyuanJiang/FitDiT) · [Leffa (MIT)](https://github.com/franciszzj/Leffa) · [identity-degradation benchmark](https://arxiv.org/pdf/2604.19748) · [imgly AGPL](https://github.com/imgly/background-removal-js/blob/main/LICENSE.md) · [RMBG-1.4 license](https://huggingface.co/briaai/RMBG-1.4) · [BiRefNet MIT](https://github.com/ZhengPeng7/BiRefNet) · [SMPL-X license](https://smpl-x.is.tue.mpg.de/modellicense.html) · [RPM shutdown](https://variety.com/2025/digital/news/netflix-acquires-ready-player-me-games-avatar-creation-1236612915/)

Commerce: [Amazon OA](https://affiliate-program.amazon.com/help/operating/agreement) · [Creators API transition](https://blog.freshstore.com/amazon-creators-api-pa-api-retirement/) · [ShareASale→Awin](https://www.awin.com/us/news-and-events/awin-news/shareasale-to-awin-upgrade) · [Awin product feeds](https://help.awin.com/developers/docs/product-feed-publisher-guide-intro) · [Impact catalogs](https://help.impact.com/partner/what-would-you-like-to-learn-about/platform-features/marketing-content/product-marketplace-and-catalogs/download-product-catalogs-as-a-partner) · [Sovrn approval](https://www.sovrn.com/blog/getting-approved-for-commerce/) · [Sovrn product API](https://developer.sovrn.com/reference/product-affiliate-api) · [Rakuten publisher reqs](https://pubhelp.rakutenadvertising.com/hc/en-us/articles/360060314292-Requirements-for-New-Publishers) · [affiliate benchmarks](https://www.dollarpocket.com/affiliate-marketing-benchmarks-2026)

Legal: [BIPA SB 2979 analysis](https://www.gtlaw.com/en/insights/2024/8/bipa-update-illinois-limits-liability-and-clarifies-electronic-consent-for-biometric-data-collection) · [try-on retailer suits](https://www.afslaw.com/perspectives/alerts/biometric-privacy-class-actions-take-aim-virtual-try-retailers) · [Lensa BIPA suit](https://www.thefashionlaw.com/ai-avatar-generator-lensa-lands-on-receiving-end-of-biometric-data-lawsuit/) · [Texas Meta settlement](https://www.texasattorneygeneral.gov/news/releases/attorney-general-ken-paxton-secures-14-billion-settlement-meta-over-its-unauthorized-capture) · [Google TX settlement](https://www.bracewell.com/resources/google-agrees-to-1-375-billion-settlement-as-texas-attorney-general-continues-data-privacy-push/) · [WMHMD biometrics](https://hintzelaw.com/blog/wa-my-health-my-data-act-pt7-biometrics) · [CO HB 24-1130](https://www.insideprivacy.com/data-privacy/colorado-privacy-act-amended-to-include-biometric-data-provisions/) · [GDPR Art. 9](https://gdpr-info.eu/art-9-gdpr/) · [COPPA 2025 rule](https://www.federalregister.gov/documents/2025/04/22/2025-05904/childrens-online-privacy-protection-rule) · [nominative fair use](https://www.inta.org/fact-sheets/fair-use-of-trademarks-intended-for-a-non-legal-audience/) · [LLC fees by state](https://www.llcuniversity.com/llc-filing-fees-by-state/) · [startup insurance costs](https://www.vouch.us/blog/startup-insurance-costs)

Competitive: [Doji raise](https://techcrunch.com/2025/05/15/doji-raises-14m-to-make-virtual-try-ons-fun-through-ai-avatars) · [Whering $7M](https://tech.eu/2026/07/07/whering-lands-7m-as-digital-wardrobe-platform-reaches-10m-users/) · [Doppl shutdown](https://support.google.com/labs/answer/16537062) · [Google try-on expansion](https://techcrunch.com/2025/10/08/googles-virtual-try-on-shopping-tool-expands-to-more-countries-now-lets-you-try-on-shoes/) · [Prime TBYB shutdown](https://www.cnbc.com/2025/01/10/amazon-to-shut-down-stitch-fix-competitor-prime-try-before-you-buy.html) · [Walmart BYOM](https://corporate.walmart.com/news/2022/09/15/walmart-levels-up-virtual-try-on-for-apparel-with-be-your-own-model-experience) · [Alta](https://techcrunch.com/2025/06/16/alta-raises-11m-to-bring-clueless-fashion-tech-to-life-with-all-star-investors/) · [Daydream](https://techcrunch.com/2025/06/25/julie-bornsteins-daydream-is-releasing-an-ai-powered-chatbot-for-fashion-related-shopping) · [Phia](https://techcrunch.com/2026/01/27/phoebe-gates-and-sophia-kiannis-phia-raises-35m-to-make-shopping-fun-again/)
